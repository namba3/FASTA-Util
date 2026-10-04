use std::time::Duration;

pub const DEFAULT_INPUT_SIZE: usize = 10_000;
pub const DEFAULT_SAMPLE_MS: u64 = 200;
pub const INVALID_BASES: &[u8] = b"xyz0123?";

pub struct Config {
    pub input_size: usize,
    pub sample_duration: Duration,
}

pub fn usage() -> &'static str {
    "Usage: cargo bench --bench nucleic_acid -- [--input-size BYTES] [--sample-ms MS]\n\
     Defaults: --input-size 10000 --sample-ms 200"
}

pub fn parse_config<I>(args: I) -> Result<Option<Config>, String>
where
    I: IntoIterator<Item = String>,
{
    let mut input_size = DEFAULT_INPUT_SIZE;
    let mut sample_ms = DEFAULT_SAMPLE_MS;
    let mut args = args.into_iter();

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--input-size" => {
                let value = args
                    .next()
                    .ok_or_else(|| format!("missing value for {arg}\n{}", usage()))?;
                input_size = value
                    .parse()
                    .map_err(|_| format!("invalid byte count '{value}' for {arg}\n{}", usage()))?;
                if input_size == 0 {
                    return Err(format!("{arg} must be greater than zero\n{}", usage()));
                }
            }
            "--sample-ms" => {
                let value = args
                    .next()
                    .ok_or_else(|| format!("missing value for {arg}\n{}", usage()))?;
                sample_ms = value
                    .parse()
                    .map_err(|_| format!("invalid duration '{value}' for {arg}\n{}", usage()))?;
                if sample_ms == 0 {
                    return Err(format!("{arg} must be greater than zero\n{}", usage()));
                }
            }
            // Cargo appends this flag when it launches a custom benchmark target.
            "--bench" => {}
            "-h" | "--help" => return Ok(None),
            _ => return Err(format!("unknown argument '{arg}'\n{}", usage())),
        }
    }

    Ok(Some(Config {
        input_size,
        sample_duration: Duration::from_millis(sample_ms),
    }))
}

pub fn repeated_bytes(bytes: &[u8], input_size: usize) -> Vec<u8> {
    (0..input_size)
        .map(|index| bytes[index % bytes.len()])
        .collect()
}

pub fn mixed_sequence(
    valid_symbols: &[u8],
    valid_percent: u32,
    mut state: u32,
    input_size: usize,
) -> Vec<u8> {
    (0..input_size)
        .map(|_| {
            // Xorshift32 keeps the generated input deterministic without a dependency.
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;

            if state % 100 < valid_percent {
                valid_symbols[(state as usize >> 8) % valid_symbols.len()]
            } else {
                INVALID_BASES[(state as usize >> 8) % INVALID_BASES.len()]
            }
        })
        .collect()
}
