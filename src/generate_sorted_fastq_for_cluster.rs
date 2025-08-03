use crate::seeding_and_filtering_seeds;
use crate::structs::{FastqRecord_isoncl_init, Minimizer_hashed};
use crate::write_output;
use crate::write_output::path_exists;
use rayon::prelude::*;
use std::collections::VecDeque;
use std::path::Path;
use std::time::Instant;
//use crate::bio_rust_file_read;
use bio::io::fastq;
use log::{debug, info};
use minimizer_iter::MinimizerBuilder;
use rustc_hash::FxHashMap;
use std::fs;

//https://doc.rust-lang.org/std/primitive.char.html#method.decode_utf16  for parsing of quality values
fn compress_sequence(seq: &[u8]) -> String {
    //compresses the sequence seq by keeping only the first character of each consecutive group of equal characters. The resulting compressed sequence is stored in the variable seq_hpol_comp.

    let mut seq_hpol_comp = String::new();
    let mut last_char: Option<&u8> = None;
    for ch in seq {
        if last_char.is_none() || last_char.unwrap() != ch {
            seq_hpol_comp.push(*ch as char);
        }
        last_char = Some(ch);
    }

    seq_hpol_comp
}

fn expected_number_errornous_kmers(quality_string: &str, k: usize, d: &[f64; 128]) -> f64 {
    //computes the expected number of errornous kmers for a read by analysing the quality entry
    let prob_error: Vec<f64> = quality_string
        .chars()
        .map(|char_| d[char_ as u8 as usize])
        .collect();
    let mut sum_of_expectations = 0.0;
    let mut current_prob_no_error = 1.0;
    //holds k nucleotides to represent a kmer
    let mut window: VecDeque<f64> = VecDeque::with_capacity(k);
    //iterates over the quality values
    for (i, &p_e) in prob_error.iter().enumerate() {
        //the probability that we do not have an error is multiplied by the probability that we did not have an error at position p_e
        current_prob_no_error *= 1.0 - p_e;

        if i >= k {
            let p_to_leave = window.pop_front().unwrap();
            current_prob_no_error /= p_to_leave;
        }

        sum_of_expectations += current_prob_no_error;

        if i >= k - 1 {
            window.push_back(1.0 - p_e);
        }
    }
    //
    (quality_string.len() - k + 1) as f64 - sum_of_expectations
}

fn calculate_error_rate(qual: &[u8], d_no_min: &[f64; 128]) -> f64 {
    let mut counts = vec![0; 128];
    let mut total_count = 0;
    let mut poisson_mean = 0.0;

    for &char_byte in qual.iter() {
        counts[char_byte as usize] += 1;
    }

    for (idx, cnt) in counts.iter().enumerate() {
        poisson_mean += *cnt as f64 * d_no_min[idx];
        total_count += *cnt;
    }

    poisson_mean / total_count as f64
}

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

//D_no_min = {chr(i) : 10**( - (ord(chr(i)) - 33)/10.0 )  for i in range(128)}
fn compute_d_no_min() -> [f64; 128] {
    let mut d = [0.0; 128];

    for i in 0..128 {
        let chr_i = i as u8 as char;
        let ord_i = chr_i as i8;
        let exponent = -(ord_i - 33) as f64 / 10.0;
        d[i] = (10.0_f64).powf(exponent);
    }
    d
}

fn analyse_fastq_and_sort(
    k: usize,
    q_threshold: f64,
    in_file_path: &str,
    quality_threshold: &f64,
    window_size: usize,
    score_vec: &mut Vec<(i32, usize)>,
    id_map: &mut FxHashMap<i32, String>,
    seeding: &str,
    s: usize,
    t: usize,
    noncanonical_bool: bool,
    verbose: bool,
) {
    /*
    Reads, filters and sorts reads from a fastq file so that we are left with reads having a reasonable quality score, that are sorted by score
     */
    let d_no_min = compute_d_no_min();
    //read_id holds the internal id we appoint to a read
    let mut read_id = 0;
    //generate a Reader object that parses the fastq-file (taken from rust-bio)
    let reader =
        fastq::Reader::from_file(Path::new(&in_file_path)).expect("We expect the file to exist");
    //make sure that we have suitable values for k_size and w_size (w_size should be larger)
    let mut w;
    if window_size > k {
        w = window_size - k + 1; // the minimizer generator will panic if w is even to break ties
        if w % 2 == 0 {
            w += 1;
        }
    }
    //k_size was chosen larger than w_size. To not fail we use every k-mer as minimizer (maybe have an error message?)
    else {
        w = 1;
    }
    //iterate over the records
    for record in reader.records() {
        let seq_rec = record.expect("invalid record");
        let header_new = seq_rec.id();
        let sequence = seq_rec.seq();
        let quality = seq_rec.qual(); //.expect("We also should have a quality");
                                      //add the read id and the real header to id_map
        if sequence.len() >= 2 * k && compress_sequence(sequence).len() >= k {
            let mut this_minimizers = vec![];
            let mut filtered_minis = vec![];
            if seeding == "minimizer" {
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

                    debug!("minimizers NEW len: {:?}", this_minimizers.len());

                    //generate_sorted_fastq_new_version::get_canonical_kmer_minimizers_hashed(sequence, k, window_size, &mut this_minimizers);
                    debug!("minimizers OLD len: {:?}", &this_minimizers.len());
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
            } else {
                debug!("No seeding");
            }
            seeding_and_filtering_seeds::filter_seeds_by_quality(
                &this_minimizers,
                quality,
                k,
                d_no_min,
                &mut filtered_minis,
                quality_threshold,
                false,
            );
            let error_rate = calculate_error_rate(quality, &d_no_min);
            if 10.0_f64 * -error_rate.log(10.0_f64) > q_threshold {
                id_map.insert(read_id, header_new.to_string());
                let score_tuple = (read_id, filtered_minis.len());
                score_vec.push(score_tuple);
                read_id += 1;
            }
        }
    }
    //sort the score_vec by the number of high-confidence seeds (most to least)
    score_vec.par_sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap()); //TODO: replace by par_sort_by
    if verbose {
        let print_vec = &score_vec[0..5];
        for score_tup in print_vec {
            debug!("ID {} count {}", &score_tup.0, score_tup.1);
        }
    }

    info!("{} reads accepted", score_vec.len());
    debug!("{:?}", score_vec.pop());
}

fn print_statistics(fastq_records: &[FastqRecord_isoncl_init]) {
    /*
    Prints the statistics for the resulting file TODO: add median and mean
     */
    let min_e = fastq_records[0].get_err_rate();
    let max_e = fastq_records[fastq_records.len() - 1].get_err_rate();
    info!("Lowest read error rate: {}", min_e);
    info!("Highest read error rate: {}", max_e);
    //logfile.write("Median read error rate:{0}\n".format(median_e))
    //logfile.write("Mean read error rate:{0}\n".format(mean_e))
}

pub(crate) fn sort_fastq_for_cluster(
    k: usize,
    q_threshold: f64,
    in_file_path: &str,
    outfolder: &String,
    quality_threshold: &f64,
    window_size: usize,
    seeding: &str,
    s: usize,
    t: usize,
    noncanonical_bool: bool,
    verbose: bool,
) {
    info!("Sorting the fastq_file");
    let now = Instant::now();
    //holds the internal ids and scores as tuples to be able to sort properly
    let mut score_vec = vec![];
    //holds the internal read id
    let mut id_map = FxHashMap::default();
    //the main step of the sort_fastq_for_cluster step: Gets the number of high-confidence seeds for each read and writes them into score_vec
    analyse_fastq_and_sort(
        k,
        q_threshold,
        in_file_path,
        quality_threshold,
        window_size,
        &mut score_vec,
        &mut id_map,
        seeding,
        s,
        t,
        noncanonical_bool,
        verbose,
    );
    let elapsed = now.elapsed();
    info!("Elapsed: {:.2?}", elapsed);
    if !path_exists(outfolder) {
        fs::create_dir(outfolder.clone()).expect("We should be able to create the directory");
    }
    //write a fastq-file that contains the reordered reads
    write_output::write_ordered_fastq(&score_vec, outfolder, &id_map, in_file_path);
    info!("Wrote the sorted fastq file");
    //print_statistics(fastq_records.borrow());
}
