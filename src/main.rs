use rand::rngs::StdRng;
use rand::SeedableRng;

mod ark;
mod zkcrypto;

fn main() {
    let rng = &mut StdRng::seed_from_u64(2322u64);
    let t = 1 << 3;

    println!("ark {}", ark::kzg10(rng, t));
    println!("zkcrypto {}", zkcrypto::kzg10(rng, t));
}
