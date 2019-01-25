public class Experiment {

	public static void main(String args[]) {

		// Compare the lat/long values for the two ellipsoids for the RGO
		Converter cAiry = new Converter(Converter.Ellipsoid.AIRY_1830, Converter.Projection.NATIONAL_GRID);
		Converter cWGS84 = new Converter(Converter.Ellipsoid.WGS_84, Converter.Projection.NATIONAL_GRID);
		
		boolean useDiagnostics = false;
		cAiry.setDiagnostics(useDiagnostics);
		cWGS84.setDiagnostics(useDiagnostics);
		
		// Figures read of screen using OS Maps for Easting/Northing of the RGO
		// Square = TQ = 51
		double RGO_E_Basis = 500000 + 38886;
		double RGO_N = 100000 + 77332;
		
		for(int offset = -15; offset <= 5; offset++) {
				
			System.out.println();
			ENPoint pNG = new ENPoint(RGO_E_Basis + offset, RGO_N);
			
			LatLongPoint pAiry = cAiry.toLatLong(pNG);
			LatLongPoint pWGS84 = cWGS84.toLatLong(pNG);
			System.out.println();
			System.out.println("++++++++++++++++++++++++++++++++++++++++++++++");
			System.out.println();
			System.out.printf(" %s %n", pNG.asString() );
			System.out.printf(" Airy:  %s %n", pAiry.asDMSString() );
			System.out.printf(" WGS84: %s %n", pWGS84.asDMSString() );
			System.out.printf(" lat diff = %f seconds %n", (pAiry.m_latitude.asDegrees() - pWGS84.m_latitude.asDegrees()) * 60 * 60 );
			System.out.printf(" lon diff = %f seconds %n", (pAiry.m_longitude.asDegrees() - pWGS84.m_longitude.asDegrees()) * 60 * 60 );
		}
		
		/*
		  Above iterations suggest that 0 longitude line for the above fixed Northing value is at about Easting 538874 using the Airy Ellipse.
		  
		  E=538874.000 / N=177332.000 
			 Airy:  Lat = 51° 28′ 38.574636″ N  Long = 0° 00′ 0.028807″ W 
			 WGS84: Lat = 51° 28′ 37.803769″ N  Long = 0° 00′ 0.762855″ W 
			 lat diff = 0.770866 seconds 
			 lon diff = 0.734048 seconds 

		  and that using WGS84, the 0 longitude line is at about Easting 538889, so about 15 metres apart ????

		  E=538889.000 / N=177332.000 
			 Airy:  Lat = 51° 28′ 38.561376″ N  Long = 0° 00′ 0.748364″ E 
			 WGS84: Lat = 51° 28′ 37.790512″ N  Long = 0° 00′ 0.014236″ E 
			 lat diff = 0.770864 seconds 
			 lon diff = 0.734127 seconds 
		*/		
	}
}
