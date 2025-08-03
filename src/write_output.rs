use crate::structs::FastqRecord_isoncl_init;
use crate::{file_actions, Cluster_ID_Map};
use log::debug;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::path::PathBuf;

pub(crate) fn write_ordered_fastq(
    score_vec: &[(i32, usize)],
    outfolder: &String,
    id_map: &FxHashMap<i32, String>,
    fastq: &str,
) {
    //writes the fastq file
    let _ = fs::create_dir_all(PathBuf::from(outfolder).join("clustering"));
    let fastq_file = File::open(fastq).unwrap();
    let mut fastq_records = FxHashMap::default();
    file_actions::parse_fastq_hashmap(fastq_file, &mut fastq_records);
    let f = File::create(outfolder.to_owned() + "/clustering/sorted.fastq")
        .expect("Unable to create file");
    let mut buf_write = BufWriter::new(&f);

    //for record in fastq_records {
    for score_tup in score_vec.iter() {
        let this_header = id_map.get(&score_tup.0).unwrap();
        let record = fastq_records.get(this_header).unwrap();
        if &record.header == this_header {
            write!(
                buf_write,
                "@{}\n{}\n+\n{}\n",
                record.get_header(),
                record.get_sequence(),
                record.get_quality()
            )
            .expect("Could not write file");
        }
    }
    buf_write.flush().expect("Failed to flush the buffer");
}

fn write_final_clusters_tsv(
    outfolder: &Path,
    clusters: Cluster_ID_Map,
    id_map: FxHashMap<i32, String>,
    header_cluster_map: &mut FxHashMap<String, i32>,
) {
    let file_path = PathBuf::from(outfolder).join("final_clusters.tsv");

    let f = File::create(file_path).expect("unable to create file");
    let mut buf_write = BufWriter::new(&f);
    let mut nr_reads = 0;
    println!("{} different clusters identified", clusters.len());
    //let nr_clusters=clusters.len();
    for (cl_id, r_int_ids) in clusters.into_iter() {
        debug!("cl_id {}, nr_reads {:?}", cl_id, nr_reads);
        for r_int_id in r_int_ids {
            let read_id = id_map.get(&r_int_id).unwrap();
            nr_reads += 1;
            let _ = writeln!(buf_write, "{}\t{}", cl_id, read_id);
            header_cluster_map.insert(read_id.to_string(), cl_id);
        }
    }
    // Flush the buffer to ensure all data is written to the underlying file
    buf_write.flush().expect("Failed to flush the buffer");
    // debug!("{} different clusters identified", nr_clusters);
    debug!("HCLM {:?}", header_cluster_map);
    println!("{} reads added to tsv file", nr_reads);
}

//TODO: this is the current RAM bottleneck: we read the whole file to then have the reads when we write the output
//Outline: sort the fastq file by cluster and then write the entries from the sorted fastq file to not having to read the full file
fn create_final_ds(
    header_cluster_map: FxHashMap<String, i32>,
    fastq: String,
    cluster_map: &mut FxHashMap<i32, Vec<FastqRecord_isoncl_init>>,
) {
    let fastq_file = File::open(fastq).unwrap();
    let mut fastq_vec = vec![];
    //parse the fastq file to store the data in fastq_vec
    file_actions::parse_fastq(fastq_file, &mut fastq_vec);
    //iterate over fastq_vec and add the reads to cluster_map
    for read in fastq_vec {
        let id = read.header.clone();
        if header_cluster_map.contains_key(&id) {
            let cluster_id = header_cluster_map.get(&id).unwrap();
            if cluster_map.contains_key(cluster_id) {
                let mut id_vec: &mut Vec<FastqRecord_isoncl_init> =
                    cluster_map.get_mut(cluster_id).unwrap();
                id_vec.push(read)
            } else {
                let id_vec = vec![read];
                cluster_map.insert(*cluster_id, id_vec);
            }
        }
    }
}

fn write_fastq_files(
    outfolder: &Path,
    cluster_map: FxHashMap<i32, Vec<FastqRecord_isoncl_init>>,
    n: usize,
) {
    let mut new_cl_id = 0;
    let mut read_cter = 0;
    //fs::create_dir_all(PathBuf::from(outfolder).join("fastq_files"));
    let fastq_outfolder = PathBuf::from(outfolder);
    //Writes the fastq files using the data structure cluster_map HashMap<i32, Vec<FastqRecord_isoncl_init>>
    for (cl_id, records) in cluster_map.into_iter() {
        if records.len() >= n {
            //only write the records if we have n or more reads supporting the cluster
            let filename = new_cl_id.to_string() + ".fastq";
            let file_path = fastq_outfolder.join(filename);
            let f = File::create(file_path).expect("unable to create file");
            let mut buf_write = BufWriter::new(&f);
            for record in records {
                write!(
                    buf_write,
                    "@{}\n{}\n+\n{}\n",
                    record.header, record.sequence, record.quality
                )
                .expect("We should be able to write the entries");
                read_cter += 1;
            }
            buf_write.flush().expect("Failed to flush the buffer");
            new_cl_id += 1; //this is the new cl_id as we skip some on the way
        }

        debug!("cl id for writing: {}, {}", cl_id, read_cter);
    }
    println!("{} reads written", read_cter);
}

pub fn path_exists(path: &str) -> bool {
    fs::metadata(path).is_ok()
}

pub(crate) fn write_output(
    outfolder: String,
    clusters: &Cluster_ID_Map,
    fastq: String,
    id_map: FxHashMap<i32, String>,
    n: usize,
    no_fastq: bool,
) {
    if !path_exists(&outfolder) {
        fs::create_dir(outfolder.clone()).expect("We should be able to create the directory");
    }
    //clustering_path: the outfolder of isONclust3
    let clustering_path = Path::new(&outfolder).join("clustering");
    if !clustering_path.exists() {
        let _ = fs::create_dir(clustering_path.clone());
    }
    //the subfolder of clustering in which we write the fastq_files ( as done in isONclust1)
    let fastq_path = clustering_path.join("fastq_files");
    if !fastq_path.exists() {
        let _ = fs::create_dir(fastq_path.clone());
    }
    let mut cluster_hashmap_fastq_record = FxHashMap::default();
    let mut header_cluster_map = FxHashMap::default();
    write_final_clusters_tsv(
        &clustering_path,
        clusters.clone(),
        id_map.clone(),
        &mut header_cluster_map,
    );
    //no_fastq: true -> we do not want to write the fastq files
    if !no_fastq {
        //create a data structure that we use to generate the proper fastq files
        create_final_ds(header_cluster_map, fastq, &mut cluster_hashmap_fastq_record);
        debug!("Cluster_hashmap: {}", cluster_hashmap_fastq_record.len());
        println!("Writing the fastq files");
        write_fastq_files(&fastq_path, cluster_hashmap_fastq_record, n);
    }
}
