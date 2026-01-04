use crate::structs;
use crate::structs::FastaRecord;
use crate::structs::FastqRecord_isoncl_init;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub(crate) fn parse_fasta(file_path: &str) -> Result<Vec<FastaRecord>, std::io::Error> {
    let file = File::open(file_path).unwrap();
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut current_header = String::new();
    let mut current_sequence = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // New header line
            if !current_header.is_empty() {
                // Store the previous sequence
                records.push(FastaRecord {
                    header: current_header.clone(),
                    sequence: current_sequence.concat(),
                });
                current_sequence.clear();
            }
            current_header = line[1..].trim().to_string(); // Remove '>'
        } else {
            // Sequence line
            current_sequence.push(line);
        }
    }

    // Store the last sequence
    if !current_header.is_empty() {
        records.push(FastaRecord {
            header: current_header.clone(),
            sequence: current_sequence.concat(),
        });
    }

    Ok(records)
}

fn shorten_header(header: &str) -> &str {
    let header_parts: Vec<_> = header.split(' ').collect();
    header_parts[0] as _
}

pub(crate) fn parse_fastq_hashmap(
    file: File,
    records: &mut FxHashMap<String, structs::FastqRecord_isoncl_init>,
) {
    //Parses a fastq file, returns a vector of FastqRecords
    let mut reader = BufReader::new(file);
    //let mut records=vec![];
    loop {
        let mut header = String::new();
        let header_read = reader.read_line(&mut header).expect("Should be contained");
        if header_read == 0 {
            break;
        }

        header = header.trim().to_owned();
        header = (header[1..]).to_string();
        header = shorten_header(&header).parse().unwrap();
        let mut sequence = String::new();
        let sequence_read = reader
            .read_line(&mut sequence)
            .expect("Should be contained");
        if sequence_read == 0 {
            break;
        }
        sequence = sequence.trim().to_owned();

        let mut quality_header = String::new();
        let quality_header_read = reader
            .read_line(&mut quality_header)
            .expect("Should be contained");
        if quality_header_read == 0 {
            break;
        }
        quality_header = quality_header.trim().to_owned();

        let mut quality = String::new();
        let quality_read = reader.read_line(&mut quality).expect("Should be contained");
        if quality_read == 0 {
            break;
        }
        quality = quality.trim().to_owned();
        let score = 0.0_f64;
        let error_rate = 0.0_f64;
        records.insert(
            header.clone(),
            FastqRecord_isoncl_init {
                header,
                sequence,
                /* quality_header,*/ quality,
                score,
                error_rate,
            },
        );
        //records.push();
    }
}

pub(crate) fn parse_fastq(file: File, records: &mut Vec<structs::FastqRecord_isoncl_init>) {
    //Parses a fastq file, returns a vector of FastqRecords
    let mut reader = BufReader::new(file);
    //let mut records=vec![];
    loop {
        let mut header = String::new();
        let header_read = reader.read_line(&mut header).expect("Should be contained");
        if header_read == 0 {
            break;
        }

        header = header.trim().to_owned();
        header = (header[1..]).to_string();
        header = shorten_header(&header).parse().unwrap();
        let mut sequence = String::new();
        let sequence_read = reader
            .read_line(&mut sequence)
            .expect("Should be contained");
        if sequence_read == 0 {
            break;
        }
        sequence = sequence.trim().to_owned();

        let mut quality_header = String::new();
        let quality_header_read = reader
            .read_line(&mut quality_header)
            .expect("Should be contained");
        if quality_header_read == 0 {
            break;
        }
        quality_header = quality_header.trim().to_owned();

        let mut quality = String::new();
        let quality_read = reader.read_line(&mut quality).expect("Should be contained");
        if quality_read == 0 {
            break;
        }
        quality = quality.trim().to_owned();
        let score = 0.0_f64;
        let error_rate = 0.0_f64;
        records.push(FastqRecord_isoncl_init {
            header,
            sequence,
            /* quality_header,*/ quality,
            score,
            error_rate,
        });
    }
}

/*pub(crate) fn parse_fastq(file: File) -> (Vec<FastqRecord_isoncl_init>, HashMap<i32,String>) {
    let file = File::open(file)?;
    let mut id_map=HashMap::new();
    let reader = BufReader::new(file);
    let mut records=vec![];
    let mut line_cter=0;
    let mut fastq_record:FastqRecord;
    let mut header;
    let mut sequence;
    let mut quality_header;
    let mut quality;
    for line in reader.lines().advance_by() {
        if line_cter % 4==0{

            header= line.unwrap();

        }
        else if  line_cter %4 ==1{
            sequence =line.unwrap();
        }
        else if line_cter%4 ==2{
            quality_header=line.unwrap();
        }
        else {
            quality =line.unwrap();
        }
        println!("{}", line?);
        line_cter+=1;
    }

    return (records,id_map)
}
pub(crate) fn parse_fastq(file: File) -> (Vec<FastqRecord_isoncl_init>, HashMap<i32,String>) {
    //Parses a fastq file, returns a vector of FastqRecords
    let mut reader = BufReader::new(file);
    let mut records=vec![];
    let mut id_int=0;
    let mut id_map=HashMap::new();
    loop {
        let mut header = String::new();
        let header_read = reader.read_line(&mut header);
        if header_read == 0 {
            break;
        }

        header = header.trim().to_owned();
        header = (&header[1..]).to_string();
        let mut sequence = String::new();
        let sequence_read = reader.read_line(&mut sequence);
        if sequence_read == 0 {
            break;
        }
        sequence = sequence.trim().to_owned();

        let mut quality_header = String::new();
        let quality_header_read = reader.read_line(&mut quality_header);
        if quality_header_read == 0 {
            break;
        }
        quality_header = quality_header.trim().to_owned();

        let mut quality = String::new();
        let quality_read = reader.read_line(&mut quality);
        if quality_read == 0 {
            break;
        }
        quality = quality.trim().to_owned();
        let score=0.0_f64;
        let error_rate=0.0_f64;
        let internal_id=id_int;
        id_map.insert(id_int, header);
        records.push(FastqRecord_isoncl_init { header,internal_id, sequence,/* quality_header,*/ quality, score,error_rate});
        id_int += 1;
    }
    return(records, id_map)
}*/
