use fasta_util::nucleic_acid::NUCLEIC_ACID_SET;
use rand::Rng;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().collect::<Vec<String>>();

    let size = args
        .get(1)
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(10000)
        .max(1);

    let mut rng = rand::thread_rng();
    let mut v = Vec::with_capacity(size);
    let set = &NUCLEIC_ACID_SET[..16];
    for _ in 0..size {
        let idx = rng.gen::<usize>() % set.len();
        v.push(set[idx]);
    }

    println!(">TestData 10000 random data");
    for line in v.chunks(50) {
        let line = String::from_utf8_lossy(line);
        println!("{line}");
    }

    Ok(())
}
