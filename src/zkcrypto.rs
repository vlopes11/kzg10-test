use bls12_381::{G1Affine, G1Projective, G2Affine, G2Projective, Scalar};
use ff::Field;
use group::Group;
use rand::rngs::StdRng;

pub fn kzg10(mut rng: &mut StdRng, t: u64) -> bool {
    let g = G1Projective::random(&mut rng);
    let h = G1Projective::random(&mut rng);
    let k = G2Projective::random(&mut rng);

    let alpha = Scalar::random(&mut rng);

    let k_alpha = k * alpha;

    let pk: Vec<G1Projective> = (0..=t)
        .map(|j| alpha.pow_vartime(&[j, 0, 0, 0]))
        .map(|a| g * a)
        .collect();

    let pk_p: Vec<G1Projective> = (0..=t)
        .map(|j| alpha.pow_vartime(&[j, 0, 0, 0]))
        .map(|a| h * a)
        .collect();

    // Commit
    let f: Vec<Scalar> = (0..=t).map(|_| Scalar::random(&mut rng)).collect();
    let f_p: Vec<Scalar> = (0..=t).map(|_| Scalar::random(&mut rng)).collect();

    let c_a: G1Projective = pk.iter().zip(f.iter()).map(|(pk, fj)| pk * fj).sum();
    let c_b: G1Projective = pk_p
        .iter()
        .zip(f_p.iter())
        .map(|(pk_p, fj)| pk_p * fj)
        .sum();

    let c = c_a + c_b;

    let (f_alpha, f_p_alpha) = {
        // Sanity check
        let f_alpha: Scalar = f
            .iter()
            .enumerate()
            .map(|(j, coeff)| alpha.pow_vartime(&[j as u64, 0, 0, 0]) * coeff)
            .sum();

        let f_p_alpha: Scalar = f_p
            .iter()
            .enumerate()
            .map(|(j, coeff)| alpha.pow(&[j as u64, 0, 0, 0]) * coeff)
            .sum();

        let c_q_a = g * f_alpha;
        let c_q_b = h * f_p_alpha;

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
    let i = Scalar::random(&mut rng);

    let k_i = k * i;

    let f_i: Scalar = f
        .iter()
        .enumerate()
        .map(|(j, coeff)| i.pow_vartime(&[j as u64, 0, 0, 0]) * coeff)
        .sum();

    let f_p_i: Scalar = f_p
        .iter()
        .enumerate()
        .map(|(j, coeff)| i.pow_vartime(&[j as u64, 0, 0, 0]) * coeff)
        .sum();

    // FIXME alpha shouldn't leak from setup
    let psi_i = (f_alpha - f_i) * (alpha - i).invert().unwrap();
    let psi_p_i = (f_p_alpha - f_p_i) * (alpha - i).invert().unwrap();

    let w_i_g = g * psi_i;
    let w_i_h = h * psi_p_i;

    let w_i = w_i_g + w_i_h;

    // VerifyEval
    let p = k_alpha - k_i;

    let q_a = g * f_i;
    let q_b = h * f_p_i;

    let q = q_a + q_b;

    let c = G1Affine::from(c);
    let w_i = G1Affine::from(w_i);
    let q = G1Affine::from(q);

    let k = G2Affine::from(k);
    let p = G2Affine::from(p);

    let a = bls12_381::pairing(&c, &k);
    let b = bls12_381::pairing(&w_i, &p);
    let c = bls12_381::pairing(&q, &k);

    if a != b + c {
        print!("pairing ");
        return false;
    }

    true
}
