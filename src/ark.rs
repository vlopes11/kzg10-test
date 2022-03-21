use ark_bls12_381::{Bls12_381, Fr};
use ark_bls12_381::{G1Projective, G2Projective};
use ark_ec::PairingEngine;
use ark_ff::Field;
use ark_std::UniformRand;
use rand::rngs::StdRng;

pub fn kzg10(rng: &mut StdRng, t: u64) -> bool {
    let g = G1Projective::rand(rng);
    let h = G1Projective::rand(rng);
    let k = G2Projective::rand(rng);

    let alpha = Fr::rand(rng);

    let mut k_alpha = k.clone();
    k_alpha *= alpha;

    let pk: Vec<G1Projective> = (0..=t)
        .map(|j| alpha.pow(&[j]))
        .map(|a| {
            let mut r = g.clone();
            r *= a;
            r
        })
        .collect();

    let pk_p: Vec<G1Projective> = (0..=t)
        .map(|j| alpha.pow(&[j]))
        .map(|a| {
            let mut r = h.clone();
            r *= a;
            r
        })
        .collect();

    // Commit
    let f: Vec<Fr> = (0..=t).map(|_| Fr::rand(rng)).collect();
    let f_p: Vec<Fr> = (0..=t).map(|_| Fr::rand(rng)).collect();

    let c_a: G1Projective = pk
        .iter()
        .copied()
        .zip(f.iter().copied())
        .map(|(mut pk, fj)| {
            pk *= fj;
            pk
        })
        .sum();

    let c_b: G1Projective = pk_p
        .iter()
        .copied()
        .zip(f_p.iter().copied())
        .map(|(mut pk_p, fj)| {
            pk_p *= fj;
            pk_p
        })
        .sum();

    let c = c_a + c_b;

    let (f_alpha, f_p_alpha) = {
        // Sanity check
        let f_alpha: Fr = f
            .iter()
            .enumerate()
            .map(|(j, coeff)| alpha.pow(&[j as u64]) * coeff)
            .sum();

        let f_p_alpha: Fr = f_p
            .iter()
            .enumerate()
            .map(|(j, coeff)| alpha.pow(&[j as u64]) * coeff)
            .sum();

        let mut c_q_a = g.clone();
        let mut c_q_b = h.clone();

        c_q_a *= f_alpha;
        c_q_b *= f_p_alpha;

        let c_q = c_q_a + c_q_b;

        if c != c_q {
            print!("commit ");
            return false;
        }

        (f_alpha, f_p_alpha)
    };

    // Open
    // VerifyPoly

    // CreateWitness
    let i = Fr::rand(rng);

    let mut k_i = k.clone();
    k_i *= i;

    let f_i: Fr = f
        .iter()
        .enumerate()
        .map(|(j, coeff)| i.pow(&[j as u64]) * coeff)
        .sum();

    let f_p_i: Fr = f_p
        .iter()
        .enumerate()
        .map(|(j, coeff)| i.pow(&[j as u64]) * coeff)
        .sum();

    // FIXME alpha shouldn't leak from setup
    let psi_i = (f_alpha - f_i) / (alpha - i);
    let psi_p_i = (f_p_alpha - f_p_i) / (alpha - i);

    let mut w_i_g = g.clone();
    let mut w_i_h = h.clone();

    w_i_g *= psi_i;
    w_i_h *= psi_p_i;

    let w_i = w_i_g + w_i_h;

    // VerifyEval
    let p = k_alpha - k_i;

    let mut q_a = g.clone();
    q_a *= f_i;

    let mut q_b = h.clone();
    q_b *= f_p_i;

    let q = q_a + q_b;

    let a = Bls12_381::pairing(c, k);
    let b = Bls12_381::pairing(w_i, p);
    let c = Bls12_381::pairing(q, k);

    if a != b + c {
        print!("pairing ");
        return false;
    }

    true
}
