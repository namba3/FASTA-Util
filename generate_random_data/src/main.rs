use fasta_util::nucleic_acid::NUCLEIC_ACID_SET;
use rand::{rngs::StdRng, RngExt, SeedableRng};
use std::io::{self, BufWriter, Write};

const DEFAULT_SIZE: usize = 10_000;
const USAGE: &str = "Usage: generate_random_data [SIZE] [--seed SEED]\n\
Generate a FASTA file containing a random nucleotide sequence.\n\
\n\
Arguments:\n\
  SIZE       Sequence length in bases (default: 10000; minimum: 1)\n\
  --seed     Use a reproducible random seed\n\
  -h, --help Print this help message";

fn parse_size_args<I>(args: I) -> Result<usize, io::Error>
where
    I: IntoIterator<Item = String>,
{
    let mut args = args.into_iter();
    let size = args
        .next()
        .map(|value| {
            value.parse::<usize>().map_err(|error| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid sequence size '{value}': {error}"),
                )
            })
        })
        .transpose()?
        .unwrap_or(DEFAULT_SIZE);

    if let Some(extra) = args.next() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("unexpected argument '{extra}'"),
        ));
    }

    Ok(size.max(1))
}

fn parse_generator_args<I>(args: I) -> Result<(usize, Option<u64>), io::Error>
where
    I: IntoIterator<Item = String>,
{
    let mut size_args = Vec::new();
    let mut seed = None;
    let mut args = args.into_iter();

    while let Some(arg) = args.next() {
        let seed_value = if arg == "--seed" {
            Some(args.next().ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput, "missing value for --seed")
            })?)
        } else {
            arg.strip_prefix("--seed=").map(str::to_owned)
        };

        if let Some(value) = seed_value {
            if seed.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--seed may only be specified once",
                ));
            }
            seed = Some(value.parse::<u64>().map_err(|error| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid seed '{value}': {error}"),
                )
            })?);
        } else {
            size_args.push(arg);
        }
    }

    Ok((parse_size_args(size_args)?, seed))
}

fn write_fasta<W, R>(output: &mut W, rng: &mut R, size: usize) -> io::Result<()>
where
    W: Write,
    R: RngExt,
{
    let set = &NUCLEIC_ACID_SET[..16];
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

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.iter().any(|arg| arg == "-h" || arg == "--help") {
        println!("{USAGE}");
        return Ok(());
    }
    let (size, seed) = parse_generator_args(args)?;

    let stdout = io::stdout();
    let mut output = BufWriter::new(stdout.lock());

    if let Some(seed) = seed {
        let mut rng = StdRng::seed_from_u64(seed);
        write_fasta(&mut output, &mut rng, size)?;
    } else {
        let mut rng = rand::rng();
        write_fasta(&mut output, &mut rng, size)?;
    }

    output.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{parse_size_args, DEFAULT_SIZE};

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_owned()).collect()
    }

    #[test]
    fn size_defaults_when_no_argument_is_given() {
        assert_eq!(parse_size_args(args(&[])).unwrap(), DEFAULT_SIZE);
    }

    #[test]
    fn size_uses_the_single_positional_argument() {
        assert_eq!(parse_size_args(args(&["123"])).unwrap(), 123);
    }

    #[test]
    fn zero_size_keeps_the_existing_minimum_of_one() {
        assert_eq!(parse_size_args(args(&["0"])).unwrap(), 1);
    }

    #[test]
    fn size_rejects_invalid_and_extra_arguments() {
        assert!(parse_size_args(args(&["abc"])).is_err());
        assert!(parse_size_args(args(&["123", "extra"])).is_err());
    }
}
