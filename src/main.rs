//#![allow(warnings)]
extern crate clap;
extern crate nohash_hasher;
extern crate rayon;

//use crate::generate_sorted_fastq_new_version::{filter_minimizers_by_quality, Minimizer,get_kmer_minimizers};
//use clap::{arg, command, Command};

mod clustering;
pub mod file_actions;
mod generate_sorted_fastq_for_cluster;
mod gff_handling;
mod seeding_and_filtering_seeds;
mod structs;
mod write_output;

mod Parallelization_side;

//mod isONclust;
use crate::clustering::cluster_merging;
use crate::structs::{FastqRecord_isoncl_init, Minimizer_hashed};
use std::collections::HashMap;

use clap::Parser;
use rayon::prelude::*;
use std::collections::VecDeque;
use std::time::Instant;

use std::convert::TryFrom;
use std::io::Read;
use std::path::Path;

use memory_stats::memory_stats;

use rustc_hash::{FxHashMap, FxHashSet};

use bio::io::fastq;
use log::info;
use minimizer_iter::MinimizerBuilder;
use simple_logger::SimpleLogger;

type Seed_Map = FxHashMap<u64, Vec<i32>>; // Change here to any other hash table implementation, e.g.,  HashMap<u64, Vec<i32>, nohash_hasher::BuildNoHashHasher<u64>>;
type Cluster_ID_Map = FxHashMap<i32, Vec<i32>>; //  Change here to any other hash table implementation, e.g., HashMap<i32, Vec<i32>,nohash_hasher::BuildNoHashHasher<i32>>;

fn compute_d() -> [f64; 128] {
    let mut d = [0.0; 128];

    for i in 0..128 {
        let chr_i = i as u8 as char;
        let ord_i = chr_i as i8;
        let exponent = -(ord_i - 33) as f64 / 10.0;
        d[i] = (10.0_f64).powf(exponent).min(0.79433);
    }
    d
}

fn expected_number_errornous_kmers(quality_string: &str, k: usize, d: &[f64; 128]) -> f64 {
    //computes the expeced number of errornous kmers for a read by analysing the quality entry
    let prob_error: Vec<f64> = quality_string
        .chars()
        .map(|char_| d[char_ as u8 as usize])
        .collect();
    let mut sum_of_expectations = 0.0;
    let mut qurrent_prob_no_error = 1.0;
    let mut window: VecDeque<f64> = VecDeque::with_capacity(k);

    for (i, &p_e) in prob_error.iter().enumerate() {
        qurrent_prob_no_error *= 1.0 - p_e;

        if i >= k {
            let p_to_leave = window.pop_front().unwrap();
            qurrent_prob_no_error /= p_to_leave;
        }

        sum_of_expectations += qurrent_prob_no_error;

        if i >= k - 1 {
            window.push_back(1.0 - p_e);
        }
    }
    info!("Quality string: {}", (quality_string.len() - k + 1) as f64);
    info!("SoE: {}", sum_of_expectations);
    (quality_string.len() - k + 1) as f64 - sum_of_expectations
}

fn calculate_error_rate(qual: &str, d_no_min: &[f64; 128]) -> f64 {
    let mut poisson_mean = 0.0;
    let mut total_count = 0;

    for char_ in qual.chars().collect::<FxHashSet<_>>() {
        let count = qual.chars().filter(|&c| c == char_).count();
        let index = char_ as usize;
        poisson_mean += count as f64 * d_no_min[index];
        total_count += count;
    }

    poisson_mean / total_count as f64
}

fn get_sorted_entries(
    mini_map_filtered: FxHashMap<i32, Vec<structs::Minimizer_hashed>>,
) -> Vec<(i32, Vec<structs::Minimizer_hashed>)> {
    // Sort by the length of vectors in descending order
    let mut sorted_entries: Vec<(i32, Vec<structs::Minimizer_hashed>)> =
        mini_map_filtered.into_iter().collect();
    sorted_entries.sort_by(|a, b| b.1.len().cmp(&a.1.len()).then_with(|| a.0.cmp(&b.0)));

    sorted_entries
}

fn filter_fastq_records(
    mut fastq_records: Vec<FastqRecord_isoncl_init>,
    d_no_min: [f64; 128],
    q_threshold: f64,
    k: usize,
    d: [f64; 128],
) -> Vec<FastqRecord_isoncl_init> {
    fastq_records.par_iter_mut().for_each(|fastq_record| {
        //calculate the error rate and store it in vector errors
        if fastq_record.sequence.len() > k {
            fastq_record
                .set_error_rate(calculate_error_rate(fastq_record.get_quality(), &d_no_min));
            let exp_errors_in_kmers =
                expected_number_errornous_kmers(fastq_record.get_quality(), k, &d);
            let p_no_error_in_kmers =
                1.0 - exp_errors_in_kmers / (fastq_record.get_sequence().len() - k + 1) as f64;
            //calculate the final score and add it to fastq_record (we have a dedicated field for that that was initialized with 0.0)
            fastq_record.set_score(
                p_no_error_in_kmers * ((fastq_record.get_sequence().len() - k + 1) as f64),
            )
        }
    });
    //filter out records that have a too high error rate
    fastq_records.retain(|record| 10.0_f64 * -record.get_err_rate().log(10.0_f64) > q_threshold);
    info!(
        "Nr of records that passed the filtering: {}",
        fastq_records.len()
    );
    fastq_records
}

fn convert_cl_id(v: usize) -> Option<i32> {
    i32::try_from(v).ok()
}

#[derive(Parser, Debug)]
#[command(name = "isONclust3")]
#[command(author = "Alexander J. Petri <alexjpetri@gmail.com>")]
#[command(version = "0.0.2")]
#[command(
    about = "Clustering of long-read sequencing data into gene families",
    long_about = "isONclust is a tool for clustering either PacBio Iso-Seq reads, or Oxford Nanopore reads into clusters, where each cluster represents all reads that came from a gene."
)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[arg(long, short, help = "Path to consensus fastq file(s)")]
    fastq: String,
    #[arg(
        long,
        short,
        help = "Path to initial clusters (stored in fasta format), which is required when --gff is set"
    )]
    init_cl: Option<String>,
    #[arg(short, help = "Kmer length")]
    k: Option<usize>,
    #[arg(short, help = " window size")]
    w: Option<usize>,
    #[arg(short, help = " syncmer length")]
    s: Option<usize>,
    #[arg(short, help = " minimum syncmer position")]
    t: Option<usize>,
    #[arg(long, short, help = "Path to outfolder")]
    outfolder: String,
    #[arg(
        long,
        short,
        default_value_t = 1,
        help = "Minimum number of reads for cluster"
    )]
    n: usize,
    #[arg(
        long,
        short,
        help = "Path to gff3 file (optional parameter), requires a reference added by calling --init-cl <REFERENCE.fasta>"
    )]
    gff: Option<String>,
    #[arg(long, help = "we do not want to use canonical seeds")]
    noncanonical: bool,
    #[arg(long, help = "Run mode of isONclust (pacbio or ont")]
    mode: String,
    //TODO:add argument telling us whether we want to use syncmers instead of kmers, maybe also add argument determining whether we want to use canonical_minimizers
    #[arg(long, help = "seeding approach we choose")]
    seeding: Option<String>,
    #[arg(long, help = "quality threshold used for the data (standard: 0.9) ")]
    quality_threshold: Option<f64>,
    #[arg(long, help = "print additional information")]
    verbose: bool,
    #[arg(
        long,
        help = "Run the post clustering step during the analysis (small improvement for results but much higher runtime)"
    )]
    post_cluster: bool,
    #[arg(
        long,
        help = "Do not write the fastq_files (no write_fastq in isONclust1)"
    )]
    no_fastq: bool,
    #[arg(
        long,
        help = "Minimum overlap threshold for reads to be clustered together (Experimental parameter)"
    )]
    min_shared_minis: Option<f64>,
    #[arg(
        long,
        help = "Minimum thresholds of shared HCS for clusters to be merged during cluster merging"
    )]
    cm_mini: Option<f64>,
}

fn main() {
    //#################################################################################################
    //INITIALIZATION
    //#################################################################################################

    let cli = Cli::parse();

    // TODO - move later on.
    let level = match cli.verbose {
        true => log::LevelFilter::Debug,
        false => log::LevelFilter::Info,
    };

    SimpleLogger::new()
        .with_level(level)
        .init()
        .expect("Failed to initialize logger.");

    info!("n: {:?}", cli.n);
    info!("outfolder {:?}", cli.outfolder);

    let mode = cli.mode;
    let n = cli.n;
    let mut k;
    let mut w;
    let mut s;
    let mut t;
    let mut quality_threshold;
    let mut min_shared_minis;
    let mut cm_mini;
    //right now we only have two modes( custom settings for variables k, w, s, and t: 'ont' for reads with  3% error rate or more and 'pacbio' for reads with less than 3% error rate)
    if mode == "ont" {
        k = 13;
        w = 21;
        quality_threshold = 0.9_f64.powi(k as i32); //TODO: standard: 0.9_f64
        min_shared_minis = 0.5;
        cm_mini = 0.5;
        s = 9;
        t = 2;
    } else if mode == "pacbio" {
        k = 15;
        w = 51;
        quality_threshold = 0.98_f64.powi(k as i32); //TODO://standard: 0.97_f64
        min_shared_minis = 0.5;
        cm_mini = 0.8;
        s = 9;
        t = 3;
    } else {
        if cli.cm_mini.is_some() {
            let cm_mini = cli.cm_mini;
        }
        if cli.quality_threshold.is_some() {
            let qt = cli.quality_threshold.unwrap();
            if cli.k.is_some() {
                k = cli.k.unwrap();
            } else {
                panic!("Please set k")
            }
            w = 0;
            t = 0;
            s = 0;
            quality_threshold = qt.powi(k as i32);
            min_shared_minis = 0.5;
            cm_mini = 0.5;
        } else {
            panic!("Please set the quality_threshold")
        }
    }
    let verbose = cli.verbose;
    /*let mut verbose = false;
    if let Some(verb) = verbo{
        verbose = true;
    }*/
    if cli.quality_threshold.is_some() {
        let qt = cli.quality_threshold.unwrap();
        quality_threshold = qt.powi(k as i32);
    }
    let post_cluster = cli.post_cluster;
    /*let mut post_cluster = false;
    if let Some(npc) = no_pc{
        post_cluster = true;
    }*/

    let no_fastq = cli.no_fastq;
    /*let mut no_fastq = false;
    if let Some(nfq) = no_fq{
        no_fastq = true;
    }*/

    let noncanonical_bool = cli.noncanonical;
    /*let mut noncanonical_bool= false;
    if let Some(noncanonical)= noncan{
        noncanonical_bool = true;
    }*/
    let seeding_input = cli.seeding.as_deref();
    let mut seeding = "minimizer";
    if let Some(seed) = seeding_input {
        seeding = seed;
    }
    if cli.k.is_some() {
        k = cli.k.unwrap();
    }

    if cli.min_shared_minis.is_some() {
        min_shared_minis = cli.min_shared_minis.unwrap()
    }
    if seeding == "syncmer" {
        if cli.s.is_some() {
            s = cli.s.unwrap();
            if (k - s + 1) % 2 != 0 {
                panic!("Please set k and s so that (k-s)+1 yields an odd number")
            }
        }
        if cli.t.is_some() {
            t = cli.t.unwrap();
            if (k - s) / 2 != t {
                panic!("Please set k,s and t to fulfill (k-s)/2=t")
            }
        }
    } else if seeding == "minimizer" {
        if cli.w.is_some() {
            w = cli.w.unwrap();
            if w < k {
                panic!("Please set w greater than k")
            } else if w % 2 == 0 {
                panic!("Please choose w to be odd")
            }
        }
    }

    info!("k: {}, w: {}, s: {}, t: {}", k, w, s, t);
    info!("quality_threshold {:?}", quality_threshold);
    info!("Min shared minis: {}", min_shared_minis);

    //let k = cli.k;
    let outfolder = cli.outfolder;
    let gff_path = cli.gff.as_deref();
    //makes the read  identifiable and gives us the possibility to only use ids during the clustering step
    let mut id_map = FxHashMap::default();
    let mut clusters: Cluster_ID_Map = HashMap::default(); //FxHashMap<i32, Vec<i32>> = FxHashMap::default();

    let mut annotation_based = false;
    let filename = outfolder.clone() + "/clustering/sorted.fastq";

    // info or debug?
    info!("Using {}s as seeds", seeding);

    let now1 = Instant::now();
    {
        //main scope (holds all the data structures that we can delete when the clustering is done
        //holds the mapping of which minimizer belongs to what clusters
        //let mut shared_seed_info: FxHashMap<i32,i32>=FxHashMap::default();
        let mut cluster_map: Seed_Map = HashMap::default(); //FxHashMap<u64, Vec<i32>> = FxHashMap::default();
        let initial_clustering_path = cli.init_cl.as_deref();
        if gff_path.is_some() {
            gff_handling::gff_based_clustering(
                gff_path,
                initial_clustering_path,
                &mut clusters,
                &mut cluster_map,
                k,
                w,
                seeding,
                s,
                t,
                noncanonical_bool,
            );
            info!(
                "{} s used for parsing the annotation information",
                now1.elapsed().as_secs()
            );
            info!("{:?}", clusters);
            annotation_based = true;
        }

        if verbose {
            if let Some(usage) = memory_stats() {
                info!("Current physical memory usage: {}", usage.physical_mem);
                info!("Current virtual memory usage: {}", usage.virtual_mem);
            } else {
                info!("Couldn't get the current memory usage :(");
            }
        }

        //#################################################################################################
        //GENERATION OF ANNOTATION BASED CLUSTERS
        //#################################################################################################

        if verbose {
            if let Some(usage) = memory_stats() {
                info!("Current physical memory usage: {}", usage.physical_mem);
                info!("Current virtual memory usage: {}", usage.virtual_mem);
            } else {
                info!("Couldn't get the current memory usage :(");
            }
        }

        //#################################################################################################
        //SORTING STEP
        //#################################################################################################

        let q_threshold = 7.0;
        //count the number of reads that were too short to be clustered
        //d_no_min contains a translation for chars into quality values
        let d_no_min = seeding_and_filtering_seeds::compute_d_no_min();
        info!("{}", filename);
        let now2 = Instant::now();
        generate_sorted_fastq_for_cluster::sort_fastq_for_cluster(
            k,
            q_threshold,
            &cli.fastq,
            &outfolder,
            &quality_threshold,
            w,
            seeding,
            s,
            t,
            noncanonical_bool,
            verbose,
        );
        let now3 = Instant::now();
        if verbose {
            info!("{} s for sorting the fastq file", now2.elapsed().as_secs());

            if let Some(usage) = memory_stats() {
                info!("Current physical memory usage: {}", usage.physical_mem);
                info!("Current virtual memory usage: {}", usage.virtual_mem);
            } else {
                info!("Couldn't get the current memory usage :(");
            }
            let now3 = Instant::now();
            info!("initial clusters {:?}", clusters);
        }
        //#################################################################################################
        //CLUSTERING STEP
        //#################################################################################################
        {
            //Clustering scope ( we define a scope so that variables die that we do not use later on)
            //the read id stores an internal id to represent our read
            let mut read_id = 0;
            //this gives the percentage of high_confidence seeds that the read has to share with a cluster to be added to it
            let reader = fastq::Reader::from_file(Path::new(&filename))
                .expect("We expect the file to exist");
            for record in reader.records() {
                let seq_rec = record.expect("invalid record");
                let header_new = seq_rec.id();
                if verbose {
                    //info!("ID: {}",header_new);
                }
                let sequence = seq_rec.seq();
                let quality = seq_rec.qual();
                //add the read id and the real header to id_map
                id_map.insert(read_id, header_new.to_string());
                let mut this_minimizers = vec![];
                let mut filtered_minis = vec![];
                if seeding == "minimizer" {
                    if w > k {
                        w = w - k + 1; // the minimizer generator will panic if w is even to break ties
                        if w % 2 == 0 {
                            w += 1;
                        }
                    }
                    if noncanonical_bool {
                        let min_iter = MinimizerBuilder::<u64, _>::new()
                            .minimizer_size(k)
                            .width((w) as u16)
                            .iter(sequence);
                        for (minimizer, position) in min_iter {
                            let mini = Minimizer_hashed {
                                sequence: minimizer,
                                position: position,
                            };
                            this_minimizers.push(mini);
                        }
                    } else {
                        let min_iter = MinimizerBuilder::<u64, _>::new()
                            .canonical()
                            .minimizer_size(k)
                            .width((w) as u16)
                            .iter(sequence);
                        for (minimizer, position, _) in min_iter {
                            let mini = Minimizer_hashed {
                                sequence: minimizer,
                                position: position,
                            };
                            this_minimizers.push(mini);
                        }
                    }
                } else if seeding == "syncmer" {
                    if noncanonical_bool {
                        seeding_and_filtering_seeds::get_kmer_syncmers(
                            sequence,
                            k,
                            s,
                            t,
                            &mut this_minimizers,
                        );
                    } else {
                        seeding_and_filtering_seeds::syncmers_canonical(
                            sequence,
                            k,
                            s,
                            t,
                            &mut this_minimizers,
                        );
                    }
                }
                seeding_and_filtering_seeds::filter_seeds_by_quality(
                    &this_minimizers,
                    quality,
                    k,
                    d_no_min,
                    &mut filtered_minis,
                    &quality_threshold,
                    verbose,
                );
                // perform the clustering step
                //info!("{} filtered_minis", filtered_minis.len());
                //info!("{} this_minimizers", this_minimizers.len());
                //info!(" ");
                let mut shared_seed_infos_norm_vec: Vec<i32> = vec![0; clusters.len()];
                clustering::cluster(
                    &filtered_minis,
                    min_shared_minis,
                    &this_minimizers,
                    &mut clusters,
                    &mut cluster_map,
                    read_id,
                    &mut shared_seed_infos_norm_vec,
                );
                read_id += 1;
                if verbose {
                    if read_id % 1000000 == 0 {
                        info!("{} reads processed", read_id);
                    }
                }
            }
            info!("Generated {} clusters from clustering", clusters.len());
            info!("Finished clustering");
            info!("{} reads used for clustering", read_id);

            if verbose {
                //info!("{} s for reading the sorted fastq file and clustering of the reads", now3.elapsed().as_secs());
            }
            if let Some(usage) = memory_stats() {
                info!("Current physical memory usage: {}", usage.physical_mem);
                info!("Current virtual memory usage: {}", usage.virtual_mem);
            } else {
                info!("Couldn't get the current memory usage :(");
            }

            //post_cluster: true -> run post_clustering
            if post_cluster {
                info!("Starting post-clustering to refine clusters");
                let now_pc = Instant::now();
                let mut shared_seed_infos_vec: Vec<i32> = vec![0; clusters.len()];

                cluster_merging(
                    &mut clusters,
                    &mut cluster_map,
                    cm_mini,
                    &mut shared_seed_infos_vec,
                    verbose,
                );
                info!("{} s for post-clustering", now_pc.elapsed().as_secs());
                info!("Got {} clusters from Post-clustering", clusters.len());
                if let Some(usage) = memory_stats() {
                    info!("Current physical memory usage: {}", usage.physical_mem);
                    info!("Current virtual memory usage: {}", usage.virtual_mem);
                } else {
                    info!("Couldn't get the current memory usage :(");
                }
            }
        }
        info!(
            "{} s for clustering and post_clustering",
            now3.elapsed().as_secs()
        );
    }

    //#################################################################################################
    //FILE OUTPUT STEP
    //#################################################################################################
    let now4 = Instant::now();
    write_output::write_output(outfolder, &clusters, filename, id_map, n, no_fastq);
    info!("{} s for file output", now4.elapsed().as_secs());
    if let Some(usage) = memory_stats() {
        info!("Current physical memory usage: {}", usage.physical_mem);
        info!("Current virtual memory usage: {}", usage.virtual_mem);
    } else {
        info!("Couldn't get the current memory usage :(");
    }
    info!("{} overall runtime", now1.elapsed().as_secs());
}
