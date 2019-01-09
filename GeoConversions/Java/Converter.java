/*
 * Convert latitude and longitude to OS National Grid Easting and Northing.
 * 
 * Using formulas from Annexe C of 'A Guide to the Coordinate Systems in Great Britain',
 * which is available at:
 * 
 * 	https://www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf
 * 
 * NB Longitude east of meridian is postive.
 * 
 * Also took some inspiration from the Javascript implementation at https://www.movable-type.co.uk/scripts/latlong-os-gridref.html, 
 * particularly use of Greek characters in variable names to match what the OS doc uses. For example
 * - φ to refer to latitude variables
 * - λ to refer to longitude variables
 * Java allows these as valid letters for use in variable naming.
 */

import static java.lang.Math.pow;
import static java.lang.Math.abs;

import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.tan;

public class Converter {

	// From Annexe A.1
	// Airy 1830 ellipsoid values, from Annexe A 
	public static double SEMI_MAJOR_AXIS = 6_377_563.396;
	public static double SEMI_MINOR_AXIS = 6_356_256.909;

	
	public enum Ellipsoid {
		AIRY_1830(6_377_563.396, 6_356_256.909);
		
		public final double a;		// SEMI_MAJOR_AXIS (m)
		public final double b;		// SEMI_MINOR_AXIS (m)
		public final double e_2;	// Eccentricty squared
		
		private Ellipsoid(double a, double b) {
			this.a = a;
			this.b = b;
			this.e_2 = ((a*a) - (b*b)) / (a*a);
		}

		public String asString() {
			return String.format("%s Ellipsoid: a=%f b=%f e_2=%f", this.name(), a, b, e_2);
		}
}

	// From Annexe A.2
	// Mercator
	public enum Projection {
		
		NATIONAL_GRID(0.9996012717, 
				PositionAngle.fromDMS(49, 0, 0, +1).radians(), PositionAngle.fromDMS(02, 0, 0, -1).radians(),
				400_000, -100_000);
		
		public final double F0;		// SCALE_FACTOR_ON_CENTRAL_MERIDIAN
		public final double φ0;		// LATITUDE_OF_TRUE_ORIGIN (Radians)
		public final double λ0; 	// LONGITUDE_OF_TRUE_ORIGIN (Radians)
		public final double E0;		// EASTING_OF_TRUE_ORIGIN
		public final double N0;		// NORTHING_OF_TRUE_ORIGIN
		
		private Projection(double F0, double φ0, double λ0, double E0, double N0) {
			this.F0 = F0;
			this.φ0 = φ0;
			this.λ0 = λ0;
			this.E0 = E0;
			this.N0 = N0;
		}

		public String asString() {
			return String.format("%s Projection: F0=%f φ0=%f λ0=%f E0=%f N0=%f", this.name(), F0, φ0, λ0, E0, N0);
		}
	}
	
	public static ENPoint toEastingsNorthings(Ellipsoid ellipsoid, Projection projection, LatLongPoint point) {
		
		// Simplify references to ellipsoid parameters
		double a = ellipsoid.a;
		//double b = ellipsoid.b;
		double e_2 = ellipsoid.e_2;
		
		// Simplify references to projection parameters
		double E0 = projection.E0;
		double N0 = projection.N0;
		double F0 = projection.F0;		
		// double φ0 = projection.φ0;
		double λ0 = projection.λ0;

		double φ = point.m_latitude.asRadians();
		double λ = point.m_longitude.asRadians();

		double sinφ = sin(φ);
		double sinφ_2 = pow(sin(φ), 2);
		double cosφ = cos(φ);
		double cosφ_3 = pow(cos(φ), 3);
		double cosφ_5 = pow(cos(φ), 5);
		double tanφ_2 = pow(tan(φ), 2);
		double tanφ_4 = pow(tan(φ), 4);
		
		double M = calculateM(ellipsoid, projection, φ);
		double ν = a * F0 * Math.pow(1 - (e_2 * sinφ_2), -0.5); 
		double ρ = a * F0 * (1 - e_2) * Math.pow(1 - (e_2 * sinφ_2), -1.5);		
		double η_2 = (ν / ρ) - 1;  
		
		double I    = M + N0;
		double II   = (ν/2)   * sinφ * cosφ;
		double III  = (ν/24)  * sinφ * cosφ_3 * (5 -     tanφ_2 + 9*η_2 );
		double IIIA = (ν/720) * sinφ * cosφ_5 * (61 - 58*tanφ_2 + tanφ_4);
		double IV   = ν       * cosφ;
		double V    = (ν/6)   * cosφ_3 * (ν/ρ -  tanφ_2);
		double VI   = (ν/120) * cosφ_5 * (5 - 18*tanφ_2 + tanφ_4 + 14*η_2 - 58*tanφ_2*η_2); 

		double λDiff = λ-λ0;
		double N = I  + II * pow(λDiff, 2) + III * pow(λDiff, 4) + IIIA * pow(λDiff, 6);
		double E = E0 + IV * pow(λDiff, 1) + V   * pow(λDiff, 3) + VI   * pow(λDiff, 5);
		
		System.out.println("Standard values:");
		System.out.println(ellipsoid.asString());
		System.out.println(projection.asString());
		System.out.println();
		System.out.println("Values derived for the specific latitude/longitude value:");
		System.out.println();
		System.out.printf("- φ    : %f %n", φ);
		System.out.printf("- λ    : %f %n", λ);
		System.out.println();
		System.out.printf("- ν    : %f %n", ν);
		System.out.printf("- ρ    : %f %n", ρ);
		System.out.printf("- ηη   : %f %n", η_2);
		System.out.printf("- M    : %f %n", M);
		System.out.printf("- I    : %f %n", I);
		System.out.printf("- II   : %f %n", II);
		System.out.printf("- III  : %f %n", III);
		System.out.printf("- IIIA : %f %n", IIIA);
		System.out.printf("- IV   : %f %n", IV);
		System.out.printf("- V    : %f %n", V);
		System.out.printf("- VI   : %f %n", VI);
		System.out.println();
		System.out.printf("- N   : %f %n", N);
		System.out.printf("- E   : %f %n", E);		
		
		ENPoint enp = new ENPoint(E, N);
		
		return enp;
	}

	public static LatLongPoint toLatLong(Ellipsoid ellipsoid, Projection projection, ENPoint point) {
		
		double M = 0;

		double a = ellipsoid.a;
		//double b = ellipsoid.b;
		double e_2 = ellipsoid.e_2;

		// Simplify references to projection parameters
		double E0 = projection.E0;
		double N0 = projection.N0;
		double F0 = projection.F0;		
		double φ0 = projection.φ0;
		double λ0 = projection.λ0;

		double E = point.m_easting;
		double N = point.m_northing;
		
		double new_φ = φ0;
		double diffCheck = 0;
		double diffLimit = 0.001 * 0.01;	// 0.01 mm
		
		int loops = 0;
		do {
			loops++;
			double old_φ = new_φ;
			new_φ = (N-N0-M)/(a*F0) + old_φ;
			M = calculateM(ellipsoid, projection, new_φ);
			
			diffCheck = abs(N-N0-M);
			System.out.printf("Loop %d: new_φ = %s M=%f diff=%f %n", loops, new_φ, M, diffCheck);
		}
		while(diffCheck >= diffLimit);

		double φ1 = new_φ;

		double sinφ = sin(φ1);
		double sinφ_2 = pow(sinφ, 2);
		double cosφ = cos(φ1);
		double tanφ = tan(φ1);
		double tanφ_2 = pow(tanφ, 2);
		double tanφ_4 = pow(tanφ, 4);
		double tanφ_6 = pow(tanφ, 6);
		
		double ν = a * F0 * Math.pow(1 - (e_2 * sinφ_2), -0.5);
		double ν_3 = Math.pow(ν,3);
		double ν_5 = Math.pow(ν,5);
		double ν_7 = Math.pow(ν,7);
		double ρ = a * F0 * (1 - e_2) * Math.pow(1 - (e_2 * sinφ_2), -1.5);		
		double η_2 = (ν / ρ) - 1;  
		
		double VII = tanφ/(2*ρ*ν);  
		double VIII = tanφ/(24*ρ*ν_3) * (5 + 3*tanφ_2 + η_2 - 9*tanφ_2*η_2);
		double IX = tanφ/(720*ρ*ν_5) * (61 + 90*tanφ_2 + 45*tanφ_4);
		
		// OS doc uses the sec function here, but as sec = 1/cos, just use 1/cos below.
		double X = 1.0/(cosφ*ν); 
		double XI = 1.0/(6*cosφ*ν_3) * (ν/ρ + 2*tanφ_2);
		double XII = 1.0/(120*cosφ*ν_5) * (5 + 28*tanφ_2 + 24*tanφ_4);
		double XIIa = 1.0/(5040*cosφ*ν_7) * (61 + 662*tanφ_2 + 1320*tanφ_4 + 720*tanφ_6);

		double Ediff = E - E0;
		double φ = φ1 - VII*pow(Ediff, 2) + VIII*pow(Ediff, 4) - IX*pow(Ediff, 6);   
		double λ = λ0 + X*Ediff - XI*pow(Ediff, 3) + XII*pow(Ediff, 5) - XIIa*pow(Ediff, 7);
				
		System.out.println();
		System.out.println("Standard values:");
		System.out.println(ellipsoid.asString());
		System.out.println(projection.asString());
		System.out.println();
		System.out.println("Values for the specific easting/norting:");
		System.out.println();
		System.out.printf("- E    : %f %n", E);
		System.out.printf("- N    : %f %n", N);
		System.out.println();
		System.out.printf("- φ1   : %e %n", φ1);
		System.out.println();
		System.out.printf("- ν    : %e %n", ν);
		System.out.printf("- ρ    : %e %n", ρ);
		System.out.printf("- ηη   : %e %n", η_2);
		System.out.printf("- VII  : %e %n", VII);
		System.out.printf("- VIII : %e %n", VIII);
		System.out.printf("- IX   : %e %n", IX);
		System.out.printf("- X    : %e %n", X);
		System.out.printf("- XI   : %e %n", XI);
		System.out.printf("- XII  : %e %n", XII);
		System.out.printf("- XIIIa: %e %n", XIIa);
		System.out.printf("- φ    : %f %n", φ);
		System.out.printf("- λ    : %f %n", λ);

		LatLong.Generator g = new LatLong.Generator(); 
		LatLong lat = g.fromRadians(LatLong.AngleType.LATITUDE, φ);
		LatLong lon = g.fromRadians(LatLong.AngleType.LONGITUDE, λ);

		LatLongPoint llp = new LatLongPoint(lat, lon);		
		return llp;
	}

	private static double calculateM(Ellipsoid ellipsoid, Projection projection, double φ) {

		// Derived ellipsoid values
		double a = ellipsoid.a;
		double b = ellipsoid.b;
		double n = (a-b)/(a+b);
		double n_2 = n*n;
		double n_3 = n*n*n;

		double F0 = projection.F0;		
		double φ0 = projection.φ0;

		double M = b * F0 * 
				(   (  ( 1 + n + (5/4)*n_2 +  (5/4)*n_3 ) * (φ-φ0)  )                             -
				    (  (   3*n +     3*n_2 + (21/8)*n_3 ) * sin(φ-φ0)     * cos(φ+φ0) )           + 
				    (  (        (15/8)*n_2 + (15/8)*n_3 ) * sin(2*(φ-φ0)) * cos(2*(φ+φ0)) )       -
				    (  (                    (35/24)*n_3 ) * sin(3*(φ-φ0)) * cos(3*(φ+φ0)) )
				);

		return M;
	}
	
	public static void main(String args[]) {

		LatLong.Generator g = new LatLong.Generator(); 
		LatLong lat = g.fromDMS(LatLong.AngleType.LATITUDE, 52, 39, 27.2531, "N");
		LatLong lon = g.fromDMS(LatLong.AngleType.LONGITUDE, 1, 43, 4.5177, "E");

		System.out.println();
		System.out.printf("Lat: %s %n" , lat.asDetailString());
		System.out.printf("Lon: %s %n" , lon.asDetailString());
		System.out.println();
		
		if(lat.isValid() && lon.isValid()) {
			LatLongPoint pIn = new LatLongPoint(lat, lon);		
			ENPoint pOut = toEastingsNorthings(Ellipsoid.AIRY_1830, Projection.NATIONAL_GRID, pIn);
			System.out.println();
			System.out.println();
			System.out.println("Converted:");
			System.out.printf(" %s %n", pIn.asString() );
			System.out.println("to:");
			System.out.printf(" %s %n", pOut.asString() );
		}
		else {
			System.out.println("Invalid point");
		}

		System.out.println();
		ENPoint p2In = new ENPoint(651409.903, 313177.270);
		
		LatLongPoint p2Out = toLatLong(Ellipsoid.AIRY_1830, Projection.NATIONAL_GRID, p2In);
		System.out.println();
		System.out.println();
		System.out.println("Converted:");
		System.out.printf("%s %n", p2In.asString() );
		System.out.println("to:");
		System.out.printf("%s %n", p2Out.asDetailString() );

	}
}

