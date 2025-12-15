use crate::structs::Minimizer_hashed;
use std::time::Instant;

use crate::{Cluster_ID_Map, Seed_Map};
use log::debug;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

pub(crate) fn reverse_complement(dna: &str) -> String {
    let reverse_complement: String = dna
        .chars()
        .rev() //TODO: test whether into_par_iter works here
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => c,
        })
        .clone()
        .collect();
    reverse_complement
}

fn calculate_shared_perc(nr_sign_minis: usize, value: i32) -> f64 {
    value as f64 / nr_sign_minis as f64
}

fn new_Fx_hashset() -> FxHashSet<i32> {
    let return_set: FxHashSet<i32> = FxHashSet::default();
    return_set
}

fn detect_whether_shared(
    min_shared_minis: f64,
    shared_seed_infos: &FxHashMap<i32, i32>,
    minimizers: &[Minimizer_hashed],
) -> (bool, i32) {
    let mut most_shared = 0.0;
    let mut shared = false;
    let mut most_shared_cluster = -1;
    let nr_minis = minimizers.len();
    let mut shared_perc: f64;
    for (key, nr_shared) in shared_seed_infos {
        //TODO: test whether into_par_iter works here
        //we have more shared minis with the cluster than our threshold and this is the cluster we share the most minimizers with
        shared_perc = calculate_shared_perc(nr_minis, *nr_shared);
        debug!(
            "shared percentage between read and cluster {} : {}",
            key, shared_perc
        );
        if shared_perc > min_shared_minis && shared_perc > most_shared {
            //} && *nr_shared >=0 {
            most_shared = shared_perc;
            most_shared_cluster = *key;
            if !shared {
                shared = true;
            }
        }
    }
    (shared, most_shared_cluster)
}

//shared_seed_infos: hashmap that holds read_id->nr shared minimizers with clusters->not updated when cluster changes!
//clustering method for the case that we do not have any annotation to compare the reads against
//shared_seed_infos: hashmap that holds read_id->nr shared minimizers with clusters->not updated when cluster changes!
//clustering method for the case that we do not have any annotation to compare the reads against
pub(crate) fn cluster(
    sign_minis: &Vec<Minimizer_hashed>,
    min_shared_minis: f64,
    minimizers: &Vec<Minimizer_hashed>,
    clusters: &mut Cluster_ID_Map,
    cluster_map: &mut Seed_Map,
    id: i32,
    shared_seed_infos_norm_vec: &mut [i32],
) {
    //clusters contains the main result we are interested in: it will contain the cluster id as key and the read_ids of reads from the cluster as value
    //cluster_map contains a hashmap in which we have a hash_value for a minimizer as key and a vector of ids as a value
    let cl_id: i32 = clusters.len() as i32;
    //we do not yet have a cluster and therefore need to fill the first read into the first
    if clusters.is_empty() {
        for minimizer in sign_minis {
            //fill cluster_map with the minimizers that we found in the first read
            cluster_map.entry(minimizer.sequence).or_default();
            let vect = cluster_map.get_mut(&minimizer.sequence).unwrap();
            if !vect.contains(&cl_id) {
                vect.push(cl_id);
            }
        }
        let id_vec = vec![id];
        clusters.insert(cl_id, id_vec);
    }
    //entry represents a read in our data
    //if we already have at least one cluster: compare the new read to the cluster(s)
    else {
        let mut most_shared_cluster = -1;
        let mut shared = false;
        debug!(
            "shared_seed_infos_norm_vec BEFORE: {:?}",
            shared_seed_infos_norm_vec
        );

        for minimizer in minimizers {
            //TODO: test whether into_par_iter works here
            //if we find the minimizer hash in cluster_map: store the clusters in belongs_to
            if let Some(belongs_to) = cluster_map.get(&minimizer.sequence) {
                //iterate over belongs_to to update the counts of shared minimizers for each cluster
                for &belong_cluster in belongs_to {
                    //TODO: test whether into_par_iter works here
                    //iterate over belongs_to to update the counts of shared minimizers for each cluster
                    shared_seed_infos_norm_vec[belong_cluster as usize] += 1;
                }
            }
        }
        // we have found clusters to compare to
        if let Some((max_cluster_id, max_shared)) = shared_seed_infos_norm_vec
            .iter()
            .enumerate()
            .max_by_key(|&(_, value)| value)
        {
            let nr_minis = minimizers.len();
            let mut shared_perc: f64;
            //we have more shared minis with the cluster than our threshold and this is the cluster we share the most minimizers with
            shared_perc = calculate_shared_perc(nr_minis, *max_shared);
            if shared_perc > min_shared_minis {
                shared = true;
                most_shared_cluster = max_cluster_id as i32;
            }
        }
        //if we have a cluster that we share enough minimizers with
        if shared {
            //add the read id to read_list
            let read_list = clusters.get_mut(&most_shared_cluster).unwrap();
            if !read_list.contains(&id) {
                read_list.push(id);
            }
            //the following for-loop updates cluster_map
            debug!(
                "PUSHING HERE most_shared_cluster: {} {}",
                most_shared_cluster, &most_shared_cluster
            );
            for sign_mini in sign_minis {
                cluster_map //TODO: test whether into_par_iter works here
                    .entry(sign_mini.sequence)
                    .or_default();
                let vect = cluster_map.get_mut(&sign_mini.sequence).unwrap();
                if !vect.contains(&most_shared_cluster) {
                    vect.push(most_shared_cluster);
                    debug!("vect: {:?}", vect);
                }
            }
        }
        //we did not find a cluster that we could put the read into-> generate a new cluster
        else {
            if !sign_minis.is_empty() {
                for sign_mini in sign_minis {
                    cluster_map //TODO: test whether into_par_iter works here
                        .entry(sign_mini.sequence)
                        .or_default();
                    let vect = cluster_map.get_mut(&sign_mini.sequence).unwrap();
                    if !vect.contains(&cl_id) {
                        debug!(" cl_id: {}  &cl_id {}", cl_id, &cl_id);
                        vect.push(cl_id);
                    }
                }
            }
            let id_vec = vec![id];
            clusters.insert(cl_id, id_vec);
        }
    }
}

//takes clusters_map as input and generates cl_set_map: a Hashmap containing the cluster id as key and a hashset of seedhashes as value.
fn generate_cluster_merging_ds(
    cl_set_map: &mut FxHashMap<i32, Vec<u64>>,
    clusters_map: &mut Seed_Map,
) {
    //cl_set_map is a hashmap with cl_id -> Hashset of seed hashes
    //iterate over clusters_map
    for (mini, vec_of_ids) in clusters_map {
        //iterate over the ids that we have stored in the value of each minimizer
        for id in vec_of_ids.iter() {
            //the cluster is already in the cluster_seeds_hash_map ->add the seed hash to the hashset, otherwise add new hashset with the seed hash
            if cl_set_map.contains_key(id) {
                cl_set_map.get_mut(id).unwrap().push(*mini);
            } else {
                let this_set: Vec<u64> = vec![*mini];
                cl_set_map.insert(*id, this_set);
            }
        }
    }
}

//helper function for the post_clustering step: Updates the 'clusters' and 'clusters_map' data structures
fn update_clusters(
    clusters: &mut Cluster_ID_Map,
    clusters_map: &mut Seed_Map,
    small_hs: &[u64],
    large_cluster_id: &i32,
    small_cluster_id: &i32,
) {
    debug!("attempt: {} into {}", small_cluster_id, large_cluster_id);
    //get the infos of clusters that belong to the two clusters we want to merge
    let small_cl_info = clusters.remove(small_cluster_id).unwrap();
    let large_cl_info = clusters.get_mut(large_cluster_id).unwrap();
    //add the reads of the small cluster into the large cluster
    large_cl_info.extend(small_cl_info);
    //clusters.remove_entry(small_cluster_id);
    //also add the hashes of the small cluster into the large cluster
    if !small_hs.is_empty() {
        for seed_hash in small_hs {
            let mut cl_vec = clusters_map.get_mut(seed_hash).unwrap();
            if !cl_vec.contains(large_cluster_id) {
                cl_vec.push(*large_cluster_id);
            }
            //delete small_cluster_id from the seed hashes, so we do not have any indication of the cluster in the ds
            cl_vec.retain(|x| *x != *small_cluster_id);
        }
    }
}

fn detect_overlaps(
    cl_set_map: &FxHashMap<i32, Vec<u64>>,
    cluster_map: &mut Seed_Map,
    merge_into: &mut Vec<(i32, i32)>,
    min_shared_minis: f64,
    small_hs: &mut FxHashSet<i32>,
    shared_seed_infos_vec: &mut [i32],
    verbose: bool,
) {
    //shared_seed_infos_vec: a vector
    let mut cl_sorted: Vec<(&i32, &Vec<u64>)> = cl_set_map.iter().collect();
    cl_sorted.sort_by_key(|&(_, v)| v.len());
    for (cl_id, hashes) in cl_sorted {
        //for (cl_id, hashes) in cl_set_map{
        //iterate over the hashes for each cl_id
        for hash in hashes.iter() {
            if let Some(belongs_to) = cluster_map.get(hash) {
                //iterate over belongs_to to update the counts of shared minimizers for each cluster
                for &belong_cluster in belongs_to {
                    //do not add an overlap with this cluster
                    if belong_cluster != *cl_id {
                        shared_seed_infos_vec[belong_cluster as usize] += 1;
                    }
                }
            }
        }
        if let Some((max_cluster_id, max_shared)) = shared_seed_infos_vec
            .iter()
            .enumerate()
            .max_by_key(|&(_, value)| value)
        {
            let nr_minis = hashes.len();
            let shared_perc: f64;
            let most_shared_cluster_id = max_cluster_id as i32;
            //calculate the percentage of shared minimizers
            shared_perc = calculate_shared_perc(nr_minis, *max_shared);
            //We only merge if we share more than min_shared_minis
            if shared_perc > min_shared_minis {
                if verbose {
                    debug!("CL_ID {}, msc {}", cl_id, most_shared_cluster_id);
                    debug!(
                        "nr_minis {}, max_shared {}, shared_perc {}",
                        nr_minis, max_shared, shared_perc
                    );
                }
                debug!("ENTERING MERGE");
                //if this cluster has less minimizers than most_shared_cluster and most_shared_cluster is not in small_hs (does not get merged into another cluster)
                if nr_minis < cl_set_map.get(&most_shared_cluster_id).unwrap().len()
                    && !small_hs.contains(&most_shared_cluster_id)
                {
                    if !merge_into.contains(&(*cl_id, most_shared_cluster_id)) {
                        //add the info to merge_into that we want to merge cl_id into most_shared_cluster
                        merge_into.push((*cl_id, most_shared_cluster_id));
                        small_hs.insert(*cl_id);
                    }
                } else {
                    //the clusters have exactly the same number of seeds
                    //if cl_id is smaller than most_shared_cluster_id (we need some kind of merging same size clusters)
                    if *cl_id < most_shared_cluster_id
                        && !small_hs.contains(&most_shared_cluster_id)
                    {
                        if !merge_into.contains(&(*cl_id, most_shared_cluster_id)) {
                            //add the info to merge_into that we want to merge cl_id into most_shared_cluster
                            merge_into.push((*cl_id, most_shared_cluster_id));
                            small_hs.insert(*cl_id);
                        }
                    }
                }
            }
        }
        // clear count vector for next cluster
        for item in shared_seed_infos_vec.iter_mut() {
            *item = 0;
        }
    }
}

fn merge_clusters_from_merge_into(
    merge_into: &mut Vec<(i32, i32)>,
    clusters_map: &mut Seed_Map,
    clusters: &mut Cluster_ID_Map,
    cl_set_map: &mut FxHashMap<i32, Vec<u64>>,
    small_hs: &FxHashSet<i32>,
    not_large: &mut FxHashSet<i32>,
) {
    debug!("Merge_into_len: {}", merge_into.len());
    for (id, value) in merge_into {
        let large_id = value;
        //we might already have deleted large_id from clusters during this iteration
        if clusters.contains_key(large_id) {
            //idea here: we merge the ids into larger clusters first, smaller clusters are still bound to merge into the new cluster later
            if !small_hs.contains(large_id) {
                //merge_clusters( clusters, clusters_map, cl_set_map, large_id, id)
                let small_hs: &Vec<u64> = cl_set_map.get(id).unwrap();
                update_clusters(clusters, clusters_map, small_hs, large_id, id);
            } else {
                not_large.insert(*large_id);
            }
        }
    }
}

pub(crate) fn cluster_merging(
    clusters: &mut Cluster_ID_Map,
    cluster_map: &mut Seed_Map,
    min_shared_minis: f64,
    shared_seed_infos_vec: &mut Vec<i32>,
    verbose: bool,
) {
    //let min_shared_minis_pc = 0.5;
    debug!("min_shared_minis_pc: {}", min_shared_minis);
    //cl_set_map is a hashmap with cl_id -> Hashset of seed hashes
    let mut cl_set_map: FxHashMap<i32, Vec<u64>> = FxHashMap::default();
    if verbose {
        debug!("Cl_set_map {:?}", cl_set_map.len());
    }
    //merge_into is a vector of a tuple(cl_id1,cl_id2)
    let mut merge_into: Vec<(i32, i32)> = vec![];
    //small_hs is a HashSet storing all cluster ids that were merged into other clusters during this iteration
    let mut small_hs: FxHashSet<i32> = FxHashSet::default();
    //used to have do-while structure
    let mut first_iter = true;
    let mut not_large = FxHashSet::default();
    //continue merging as long as we still find clusters that we may merge
    while !merge_into.is_empty() || first_iter {
        //clear merge_into as this is the indicator how often we attempt to merge further (the while loop depends on it)
        merge_into.clear();
        debug!("MI {:?}", merge_into);
        small_hs.clear();
        //set first_iter to be false to not stay in a infinity loop
        first_iter = false;
        //merge_into contains the information about which clusters to merge into which
        //generate the data structure giving us merge infos
        let now_pc1 = Instant::now();
        generate_cluster_merging_ds(&mut cl_set_map, cluster_map);
        debug!("{} s for creating the pc ds", now_pc1.elapsed().as_secs());
        debug!("Post_clustering_ds generated");
        let now_pc2 = Instant::now();
        detect_overlaps(
            &mut cl_set_map,
            cluster_map,
            &mut merge_into,
            min_shared_minis,
            &mut small_hs,
            shared_seed_infos_vec,
            verbose,
        );
        debug!(
            "{} s for detection of overlaps",
            now_pc2.elapsed().as_secs()
        );
        if verbose {
            debug!("Merge_into {:?}", merge_into);
        }
        let now_pc3 = Instant::now();
        merge_clusters_from_merge_into(
            &mut merge_into,
            cluster_map,
            clusters,
            &mut cl_set_map,
            &small_hs,
            &mut not_large,
        );
        debug!("{} s for merging of clusters", now_pc3.elapsed().as_secs());
        let now_pc4 = Instant::now();
        merge_into.retain(|&(_, second)| !not_large.contains(&second));
        debug!("{} s for retaining merge_into", now_pc4.elapsed().as_secs());
        debug!("{} s since create ds", now_pc2.elapsed().as_secs());

        cl_set_map.clear();
    }
    debug!("min_shared_minis_pc: {}", min_shared_minis);
}

pub(crate) fn generate_initial_cluster_map(
    this_minimizers: &Vec<Minimizer_hashed>,
    init_cluster_map: &mut Seed_Map,
    identifier: i32,
) {
    for minimizer in this_minimizers {
        init_cluster_map.entry(minimizer.sequence).or_default();
        // Check if id was retained (not a duplicate) and push it if needed
        let vec = init_cluster_map.get_mut(&minimizer.sequence).unwrap();
        if !vec.contains(&identifier) {
            vec.push(identifier);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_reverse_complement() {
        let rev_comp = reverse_complement("GGGGATCATCAGGGCTA");
        assert_eq!(rev_comp, "TAGCCCTGATGATCCCC");
        let rev_comp2 = reverse_complement("ATCGA");
        assert_eq!(rev_comp2, "TCGAT");
    }
}
