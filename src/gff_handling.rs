use crate::structs::{Coord_obj, Minimizer_hashed};
use std::fs::File;
use std::io::{BufRead, BufReader};

use bio::io::gff;
use rustc_hash::{FxHashMap, FxHashSet};

use bio::io::gff::GffType::GFF3;
use log::info;
use rayon::prelude::*;

use bio::io::fasta;
use minimizer_iter::MinimizerBuilder;
use std::path::Path;
extern crate rayon;
use crate::clustering;
use crate::{seeding_and_filtering_seeds, ClusterIDMap, SeedMap};
use log::debug;
use std::time::Instant;

//TODO: add overlap detection
//TODO: remove multiple occasions of minimizers in the same gene if the exons overlap
//TODO: possibly use the gene id for cluster_identification

fn detect_overlaps(
    gene_map: &FxHashMap<i32, Vec<Coord_obj>>,
    this_gene_id: &i32,
    this_coords: &Vec<Coord_obj>,
    overlap_ctr: &mut i32,
) {
    for (gene_id_other, coords) in gene_map {
        if gene_id_other > this_gene_id {
            for coord in coords {
                for this_coord in this_coords {
                    if coord.overlaps_with(this_coord) {
                        *overlap_ctr += 1;
                    }
                }
            }
        }
    }
}

fn parse_fasta_and_gen_clusters(
    fasta_path: Option<&str>,
    coords: FxHashMap<String, FxHashMap<i32, Vec<Coord_obj>>>,
    clusters: &mut ClusterIDMap,
    init_cluster_map: &mut SeedMap,
    k: usize,
    w: usize,
) {
    info!("parse_fasta");
    let path = fasta_path.unwrap();
    let reader = fasta::Reader::from_file(Path::new(path)).expect("We expect the file to exist");
    //let mut reader = parse_fastx_file(&filename).expect("valid path/file");
    reader.records().for_each(|record| {
        let mut internal_id = 0;
        let mut record_minis = vec![];
        //retreive the current record
        let seq_rec = record.expect("invalid record");
        let sequence = std::str::from_utf8(seq_rec.seq()).unwrap().to_uppercase();
        let mut overlap_ctr = 0;
        //in the next lines we make sure that we have a proper header and store it as string
        let id = seq_rec.id().to_string().split(' ').collect::<Vec<_>>()[0].to_string();
        info!("Now to the coords_ds");
        if let Some(gene_map) = coords.get(id.as_str()) {
            //iterate over the genes in the gene_map
            for (gene_id, exon_coords) in gene_map {
                //TODO:test whether into_par_iter works here
                let mut coords_in_gene = FxHashSet::default();
                for exon_coord in exon_coords {
                    if !coords_in_gene.contains(exon_coord) {
                        let exon_seq = &sequence
                            [exon_coord.startpos as usize..exon_coord.endpos as usize]
                            .to_string();
                        let mut exon_minis = vec![];
                        seeding_and_filtering_seeds::get_canonical_kmer_minimizers_hashed(
                            exon_seq.as_bytes(),
                            k,
                            w,
                            &mut exon_minis,
                        );
                        record_minis.append(&mut exon_minis);
                        coords_in_gene.insert(exon_coord);
                    }
                }
                detect_overlaps(gene_map, gene_id, exon_coords, &mut overlap_ctr);
                debug!("Record_seq {}: {:?}", id, record_minis);
                clustering::generate_initial_cluster_map(&record_minis, init_cluster_map, *gene_id);
                let id_vec = vec![];
                clusters.insert(*gene_id, id_vec);
                debug!("{:?}", init_cluster_map);
                internal_id += 1;
            }
        }
        info!("{} overlaps between genes (their exons) ", overlap_ctr);
    });
}

fn parse_gtf_and_collect_coords(
    gtf_path: Option<&str>,
    coords: &mut FxHashMap<String, FxHashMap<i32, Vec<Coord_obj>>>,
) {
    let reader = gff::Reader::from_file(gtf_path.unwrap(), GFF3);
    let mut gene_id = 0;
    let mut coords_in_gene = FxHashSet::default();
    let mut true_gene = false;
    for record in reader.expect("The reader should find records").records() {
        let rec = record.ok().expect("Error reading record.");
        //we have a new gene
        if rec.feature_type() == "gene" {
            //|| rec.feature_type() == "pseudogene"{
            true_gene = true;
            //see if we are in a new chromosome/scaffold
            if !coords.contains_key(rec.seqname()) {
                //we are in a new chromosome/scaffold
                let sname = rec.seqname().to_string();
                coords.insert(sname, FxHashMap::default());
            }
            //increase the gene_id by 1
            gene_id += 1;
            coords_in_gene.clear();
        } else if rec.feature_type() == "exon" {
            //we only are interested in exons from true genes
            if true_gene {
                debug!("{} {} {}", rec.seqname(), rec.feature_type(), gene_id);
                if let Some(gene_map) = coords.get_mut(rec.seqname()) {
                    if let Some(coord_vec) = gene_map.get_mut(&gene_id) {
                        let coord_o = Coord_obj {
                            startpos: *rec.start(),
                            endpos: *rec.end(),
                        };
                        if !coords_in_gene.contains(&coord_o) {
                            coord_vec.push(coord_o.clone());
                            coords_in_gene.insert(coord_o);
                        }
                    } else {
                        let coord_o = Coord_obj {
                            startpos: *rec.start(),
                            endpos: *rec.end(),
                        };
                        if !coords_in_gene.contains(&coord_o) {
                            let mut coord_vec = vec![];
                            coord_vec.push(coord_o.clone());
                            gene_map.insert(gene_id, coord_vec);
                            coords_in_gene.insert(coord_o);
                        }
                    }
                }
            }
        }
        //we skip any exons found in pseudogenes for now
        else if rec.feature_type() == "pseudogene" {
            true_gene = false;
        }
    }
}

pub(crate) fn resolve_gff(
    gff_path: Option<&str>,
    fasta_path: Option<&str>,
    clusters: &mut ClusterIDMap,
    cluster_map: &mut SeedMap,
    k: usize,
    w: usize,
) {
    info!("Resolving GFF file ");
    let now1 = Instant::now();
    let mut coords = FxHashMap::default(); //: HashMap<K, HashMap<i32, Vec<Coord_obj>, BuildHasherDefault<FxHasher>>, BuildHasherDefault<FxHasher>> = FxHashMap::default();
    parse_gtf_and_collect_coords(gff_path, &mut coords);
    info!(
        "{} s used for parsing the gff file",
        now1.elapsed().as_secs()
    );
    info!("First step done");
    let now2 = Instant::now();
    parse_fasta_and_gen_clusters(fasta_path, coords, clusters, cluster_map, k, w);
    info!(
        "Generated {} initial clusters from the reference",
        clusters.len()
    );
    info!(
        "{} s used for parsing the fasta file",
        now2.elapsed().as_secs()
    );
    info!("{} s for full GFF resolution", now1.elapsed().as_secs());
    info!("GTF resolved");
}

pub(crate) fn gff_based_clustering(
    gff_path: Option<&str>,
    fasta_path: Option<&str>,
    clusters: &mut ClusterIDMap,
    cluster_map: &mut SeedMap,
    k: usize,
    w: usize,
    seeding: &str,
    s: usize,
    t: usize,
    noncanonical_bool: bool,
) {
    // Read the FASTA file
    let fasta_reader = File::open(Path::new(fasta_path.unwrap())).unwrap();
    let fasta_buf_reader = BufReader::new(fasta_reader);
    let fasta_records = fasta::Reader::new(fasta_buf_reader).records();
    // Read the GFF file
    let gff_reader = File::open(Path::new(gff_path.unwrap())).unwrap();
    let gff_buf_reader = BufReader::new(gff_reader);
    let mut binding = gff::Reader::new(gff_buf_reader, GFF3);
    let mut gff_records = binding.records();
    let mut gene_id = 0;
    let mut previous_genes = 0;
    // Iterate through FASTA records
    for fasta_record in fasta_records {
        let fasta_record = fasta_record.expect("Error reading FASTA record");
        let scaffold_id = fasta_record.id().to_string();
        let sequence = std::str::from_utf8(fasta_record.seq())
            .unwrap()
            .to_uppercase();
        let record_minis = vec![];
        let mut is_gene = false;
        debug!("scaffold {}", scaffold_id);
        // Process GFF records for the current scaffold ID
        for gff_record in gff_records.by_ref() {
            let gff_record = gff_record.expect("Error reading GFF record");
            let gff_scaffold_id = gff_record.seqname().to_string();
            // Check if the scaffold IDs match
            if scaffold_id == gff_scaffold_id {
                if gff_record.feature_type() == "gene"
                    && gff_record.attributes().get("gene_biotype").expect(
                        "This algorithm requires a gene_biotype to extract the coding genes",
                    ) == "protein_coding"
                {
                    gene_id += 1;
                    is_gene = true;
                } else if gff_record.feature_type() == "exon" {
                    let exon_seq =
                        &sequence[*gff_record.start() as usize..*gff_record.end() as usize];
                    let mut this_minimizers = vec![];
                    if seeding == "minimizer" {
                        if noncanonical_bool {
                            let min_iter = MinimizerBuilder::<u64, _>::new()
                                .minimizer_size(k)
                                .width((w) as u16)
                                .iter(sequence.as_bytes());
                            for (minimizer, position) in min_iter {
                                let mini = Minimizer_hashed {
                                    sequence: minimizer,
                                    position,
                                };
                                this_minimizers.push(mini);
                            }
                        } else {
                            let min_iter = MinimizerBuilder::<u64, _>::new()
                                .canonical()
                                .minimizer_size(k)
                                .width((w) as u16)
                                .iter(sequence.as_bytes());
                            for (minimizer, position, _) in min_iter {
                                let mini = Minimizer_hashed {
                                    sequence: minimizer,
                                    position,
                                };
                                this_minimizers.push(mini);
                            }
                        }
                    } else if seeding == "syncmer" && exon_seq.len() > s {
                        seeding_and_filtering_seeds::syncmers_canonical(
                            exon_seq.as_bytes(),
                            k,
                            s,
                            t,
                            &mut this_minimizers,
                        );
                    }
                } else if gff_record.feature_type() == "pseudogene" {
                    is_gene = false;
                }
                clustering::generate_initial_cluster_map(&record_minis, cluster_map, gene_id);
                let id_vec = vec![];
                clusters.insert(gene_id, id_vec);
            } else {
                info!(
                    "found {} genes in {}",
                    gene_id - previous_genes,
                    scaffold_id
                );
                previous_genes = gene_id;
                if gff_record.feature_type() == "gene"
                    && gff_record.attributes().get("gene_biotype").expect(
                        "This algorithm requires a gene_biotype to extract the coding genes",
                    ) == "protein_coding"
                {
                    gene_id += 1;
                    is_gene = true;
                }
                // If scaffold IDs don't match, break the inner loop
                break;
            }
        }
    }
}
