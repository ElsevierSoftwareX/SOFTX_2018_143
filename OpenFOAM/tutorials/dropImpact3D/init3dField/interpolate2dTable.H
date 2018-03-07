scalar interpolate2dTable
(
	const scalar& r,
	const scalar& y,
	const Field<scalarField>& FieldData,
	const scalar& spacing,
	const scalar& defaultValue
)
{
	scalar interpolatedValue = defaultValue;

	int i = -1;
	int j = -1;

	int nDataR = FieldData.size();
	int nDataY = FieldData[0].size();

	while ((i < nDataR) && (r-SMALL >= (i+1)*spacing)){i++;}
	while ((j < nDataY) && (y-SMALL >= (j+1)*spacing)){j++;}

	if ((i <= nDataR-2) && (j <= nDataY-2) && (i >= 0) && (j >= 0))
	{
        scalar xi = r/spacing - i;
        scalar eta = y/spacing - j;

		if ((xi < -1e-3) || (eta < -1e-3) || (xi > 1.001) || (eta > 1.001))
		{
			Info << "ausserhalb des Definitionsbereichs: "  << xi << " "<<  eta << endl;
		}

		scalar w1 = (1.0-xi)*(1.0-eta);
		scalar w2 = xi*(1.0-eta);
		scalar w3 = eta*(1.0-xi);
		scalar w4 = xi*eta;

		interpolatedValue = w1*FieldData[i][j] + w2*FieldData[i+1][j] + w3*FieldData[i][j+1] + w4*FieldData[i+1][j+1];
	}
	else
	{
		if ((i <= nDataR-1) && (j <= nDataY-1) && (i >= 0) && (j >= 0))
		{
			interpolatedValue = FieldData[i][j];
		}
	}

	return interpolatedValue;
}
