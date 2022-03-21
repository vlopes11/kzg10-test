pub trait PolyCommit {
    fn setup();

    fn commit();

    fn open();

    fn verify_poly();

    fn create_witness();

    fn verify_eval();
}
