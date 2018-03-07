    // Read fields
    volScalarField alpha1
    (
        IOobject
        (
            "alpha.water",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    label nPointsR = 430;
    label nPointsZ = 1000;
    scalar resolution = 10e-6;


    //index 1 --> R
    //index 2 --> Z

    Field<scalarField> inputUrad (nPointsR,scalarField(nPointsZ,0.0));
    Field<scalarField> inputUax (nPointsR,scalarField(nPointsZ,0.0));
    Field<scalarField> inputAlpha (nPointsR,scalarField(nPointsZ,0.0));

    

	//- read data
	Info << "Reading input field for alpha1" << nl << endl;
	std::ostringstream inputStringAlpha;
	inputStringAlpha << "inputFields/alpha";
	std::ifstream inputDataAlpha;
	inputDataAlpha.open(inputStringAlpha.str().c_str());

	for (int j=0; j<nPointsZ; j++)
	{
		for (int i=0; i<nPointsR; i++)
		{
			inputDataAlpha >> inputAlpha[i][j];
		}
	}
	inputDataAlpha.close();


	Info << "Reading input field for Urad" << nl << endl;
	std::ostringstream inputStringUrad;
	inputStringUrad << "inputFields/Urad";
	std::ifstream inputDataUrad;
	inputDataUrad.open(inputStringUrad.str().c_str());

	for (int j=0; j<nPointsZ; j++)
	{
		for (int i=0; i<nPointsR; i++)
		{
			inputDataUrad >> inputUrad[i][j];
		}
	}
	inputDataUrad.close();


	Info << "Reading input field for Uax" << nl << endl;
	std::ostringstream inputStringUax;
	inputStringUax << "inputFields/Uax";
	std::ifstream inputDataUax;
	inputDataUax.open(inputStringUax.str().c_str());

	for (int j=0; j<nPointsZ; j++)
	{
		for (int i=0; i<nPointsR; i++)
		{
			inputDataUax >> inputUax[i][j];
		}
	}
	inputDataUax.close();
