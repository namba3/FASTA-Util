use fasta_util::nucleic_acid::NUCLEIC_ACID_SET;
use rand::RngExt;
use std::io::{self, BufWriter, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().collect::<Vec<String>>();

    let size = args
        .get(1)
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(10000)
        .max(1);

    let mut rng = rand::rng();
    let set = &NUCLEIC_ACID_SET[..16];

    let stdout = io::stdout();
    let mut output = BufWriter::new(stdout.lock());
    writeln!(output, ">TestData {size} random data")?;

    let mut remaining = size;
    let mut line = [0; 50];
    while remaining > 0 {
        let line_len = remaining.min(line.len());
        for base in &mut line[..line_len] {
            *base = set[rng.random_range(0..set.len())];
        }
        output.write_all(&line[..line_len])?;
        output.write_all(b"\n")?;
        remaining -= line_len;
    }

    output.flush()?;
    Ok(())
}
