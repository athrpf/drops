Directory, that contains various tests for parallel DROPS version. In
particular these tests are:

- TestExchangePar:
	o Content: Test for the ExchangeCL, which handles the communication 
		of numerical data
	o Parameter file: param-files/Exchange.param
	o Reference output for 4 processes: ref-out/Exchange_P004.txt

- TestPoissonPar:
	o Content: Test for parallel solvers applied to a discretized 
		Poisson problem
	o Parameter file: param-files/Poisson.param
	o Reference output for 4 processes: ref-out/Poisson_P004.txt

- TestRefPar:
	o Content: Test for parallel refinement-/coarsening algorithm
	o Parameter file: param-files/Ref.param
	o Reference output for 4 processes: ref-out/Ref_P004.txt

- TestInterpolPar
	o Content: Test parallel interpolation due to changed multigrid
	o Parameter file: none
	o Reference output for 4 processes: ref-out/Interpol_P004.txt

- TestSedPar:
	o Content: Test for parallel sedimenting of a droplet (just setting
		droplet to right position without performing numerics 
	o Parameter file: param-files/Sed.param
	o Reference output for 4 processes: ref-out/Sed_P004.txt

