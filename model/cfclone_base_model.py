def cfclone_base_model(ctdna: np.array, clone_cn_profiles: np.array, scale: float) -> None:
    """Inferred parameters are clone proportions.

    Args:
        ctdna: numpy array (n, ) of gc and mappability corrected ctdna binned read counts outputted from HMMcopy.
        clone_cn_profiles: numpy array (n, num_clones) of clone copy number profiles.
        scale: Scale (variance) parameter for student-t distribution
    """
    num_clones = clone_cn_profiles.shape[1]
    rho = numpyro.sample('rho', Dirichlet(jnp.ones(num_clones)))
    mu = jnp.log(jnp.sum(clone_cn_profiles*rho, axis=1)) - jnp.log(jnp.mean(jnp.sum(clone_cn_profiles*rho, axis=1)))
    with numpyro.plate('data', size=len(ctdna)):
        numpyro.sample('obs', StudentT(df=2, loc=mu, scale=scale), obs=ctdna)