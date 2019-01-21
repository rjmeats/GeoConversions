
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
 * - φ (phi) in latitude variables
 * - λ (lambda) in longitude variables
 * Java allows these as valid letters for use in variable naming.
 */

// Static imports of Math methods for concise use later.
import static java.lang.Math.pow;
import static java.lang.Math.abs;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.tan;

public class Converter {

	// From Annexe A.1
	// Airy 1830 ellipsoid values, potentially other ellipsoids
	
	public enum Ellipsoid {
		AIRY_1830(6_377_563.396, 6_356_256.909),
		WGS_84(6_378_137.000, 6_356_752.3141);
		
		public final double a;		// Semi-major axis (m)
		public final double b;		// Semi-minor axis (m)
		public final double e_2;	// Eccentricty squared
		
		private Ellipsoid(double a, double b) {
			this.a = a;
			this.b = b;
			this.e_2 = ((a*a) - (b*b)) / (a*a);
		}

		public String asString() {
			return String.format("%s Ellipsoid: a=%f b=%f e²=%f", this.name(), a, b, e_2);
		}
	}

	// ======================================================================================
	// ======================================================================================
	
	// From Annexe A.2
	// Mercator projection parameters for the National Grid, and potentially other systems
	
	public enum Projection {
		
		NATIONAL_GRID(0.9996012717, 
				PositionAngle.fromDMS(49, 0, 0, +1).radians(), PositionAngle.fromDMS(02, 0, 0, -1).radians(),
				400_000, -100_000);
		
		public final double F0;		// Scale factor on central meridian
		public final double φ0;		// Latitude of true origin (Radians)
		public final double λ0; 	// Longitude of true origin (Radians)
		public final double E0;		// Easting of true origin (m)
		public final double N0;		// Easting of true northing (m)
		
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
	
	// ======================================================================================
	// ======================================================================================

	Ellipsoid m_ellipsoid;
	Projection m_projection;

	// Control production of diagnostics to stdout
	private boolean m_diagnostics;

	// Simplified references to ellipsoid and projection values.
	private double a;
	private double b;
	private double e_2;
	
	private double E0;
	private double N0;
	private double F0;		
	private double φ0;
	private double λ0;

	Converter(Ellipsoid ellipsoid, Projection projection) {
		m_ellipsoid = ellipsoid;
		m_projection = projection;
		m_diagnostics = false;
		
		// Simplify references to ellipsoid parameters
		a = ellipsoid.a;
		b = ellipsoid.b;
		e_2 = ellipsoid.e_2;
		
		// Simplify references to projection parameters
		E0 = projection.E0;
		N0 = projection.N0;
		F0 = projection.F0;		
		φ0 = projection.φ0;
		λ0 = projection.λ0;		
	}

	void setDiagnostics(boolean onOff) {
		m_diagnostics = onOff;
	}

	// Convert from latitude/longitude to Easting/Northing for a specific ellipsoid and Mercator projection 
	public ENPoint toEastingNorthing(LatLongPoint point) {
		
		double φ = point.m_latitude.asRadians();		// phi = latitude in radians
		double λ = point.m_longitude.asRadians();		// lambda = longitude in radians 

		// Pre-calculate some trig functions and their powers of the latitude value
		double sinφ = sin(φ);
		double cosφ = cos(φ);
		double tanφ = tan(φ);
		double sinφ_2 = pow(sinφ, 2);
		double cosφ_3 = pow(cosφ, 3);
		double cosφ_5 = pow(cosφ, 5);
		double tanφ_2 = pow(tanφ, 2);
		double tanφ_4 = pow(tanφ, 4);
		
		// Some intermediate values, using variable names matching OS document. Not clear what they represent! 
		// M - something relating to the Meridian arc (northing). Complex calculation also used in converting in the other direction, so has its own method
		// ν = nu
		// ρ = rho
		// η = eta
		double M = calculateM(φ);
		double ν = a*F0 * pow(1 - (e_2 * sinφ_2), -0.5); 
		double ρ = a*F0 * (1 - e_2) * pow(1 - (e_2 * sinφ_2), -1.5);		
		double η_2 = (ν/ρ) - 1;  
		
		// Coefficients for components of the final N and E calculation. Again using OS naming.
		// Northing coefficients
		double I    = M;		// OS doc adds N0 here, but I do it below to make N and E calculations more similar.
		double II   = (ν/2)   * sinφ * cosφ;
		double III  = (ν/24)  * sinφ * cosφ_3 * (5 -     tanφ_2 + 9*η_2 );
		double IIIA = (ν/720) * sinφ * cosφ_5 * (61 - 58*tanφ_2 + tanφ_4);
		// Easting coefficients
		double IV   = ν       * cosφ;
		double V    = (ν/6)   * cosφ_3 * (ν/ρ -  tanφ_2);
		double VI   = (ν/120) * cosφ_5 * (5 - 18*tanφ_2 + tanφ_4 + 14*η_2 - 58*tanφ_2*η_2); 

		// Apply the coefficients to powers of the difference between our longitude and the projection origin longitude to obtain N and E values  
		double λDiff = λ-λ0;
		double N = N0 + I + II*pow(λDiff, 2) + III*pow(λDiff, 4) + IIIA*pow(λDiff, 6);
		double E = E0 +     IV*pow(λDiff, 1) + V*pow(λDiff, 3)   + VI*pow(λDiff, 5);
		
		if(m_diagnostics) {
			System.out.println("Standard values:");
			System.out.println(m_ellipsoid.asString());
			System.out.println(m_projection.asString());
			System.out.println();
			System.out.println("Values derived for the specific latitude/longitude value:");
			System.out.println();
			System.out.printf("- φ    : %e %n", φ);
			System.out.printf("- λ    : %e %n", λ);
			System.out.println();
			System.out.printf("- ν    : %e %n", ν);
			System.out.printf("- ρ    : %e %n", ρ);
			System.out.printf("- η²   : %e %n", η_2);
			System.out.printf("- M    : %e %n", M);
			System.out.printf("- I    : %e %n", I+N0);		// Add N0 here to match OS calculation value
			System.out.printf("- II   : %e %n", II);
			System.out.printf("- III  : %e %n", III);
			System.out.printf("- IIIA : %e %n", IIIA);
			System.out.printf("- IV   : %e %n", IV);
			System.out.printf("- V    : %e %n", V);
			System.out.printf("- VI   : %e %n", VI);
			System.out.println();
			System.out.printf("- N   : %f m%n", N);
			System.out.printf("- E   : %f m%n", E);		
		}
		
		ENPoint enp = new ENPoint(E, N);
		
		return enp;
	}

	public LatLongPoint toLatLong(ENPoint point) {
		
		double E = point.m_easting;
		double N = point.m_northing;

		// No exact formula, use an iterative approach to calculate what the latitude would be (φ1) for this
		// northing if the easting was that of the true origin (i.e. E = E0). Calculate the 'M' value for
		// the latitude until M is very close to the N-N0 northing offset.
		
		double diffLimit = 0.001 * 0.01;	// Keep iterating until difference is 0.01 mm or less		
		double M = 0;
		double φ1 = φ0;						// Start with the latitude set to that of the true origin. 
		int loops = 0;
		while(true) {
			loops++;
			double old_φ = φ1;
			φ1 = (N-N0-M)/(a*F0) + old_φ;
			M = calculateM(φ1);
			
			double diffCheck = abs(N-N0-M);
			if(m_diagnostics) {
				System.out.printf("Loop %d: new_φ = %s M=%f diff=%f %n", loops, φ1, M, diffCheck);
			}
			
			if(diffCheck < diffLimit) break;
			
			// Protect against non-convergence
			if(loops > 100) {
				System.err.println("Excessive loops calculating latitude");
				break;
			}
		}

		// Pre-calculate some trig functions and their powers of the latitude value
		double sinφ = sin(φ1);
		double cosφ = cos(φ1);
		double tanφ = tan(φ1);
		double sinφ_2 = pow(sinφ,2);
		double tanφ_2 = pow(tanφ,2);
		double tanφ_4 = pow(tanφ,4);
		double tanφ_6 = pow(tanφ,6);
		
		// As before, some intermediate values, using variable names matching OS document. Not clear what they represent! 
		// ν = nu
		// ρ = rho
		// η = eta
		double ν = a*F0 * pow(1 - (e_2 * sinφ_2), -0.5); 
		double ρ = a*F0 * (1 - e_2) * pow(1 - (e_2 * sinφ_2), -1.5);		
		double η_2 = (ν/ρ) - 1;  
		double ν_3 = pow(ν,3);
		double ν_5 = pow(ν,5);
		double ν_7 = pow(ν,7);
		
		// Coefficients for components of the final lat and long calculation. Again using OS naming.
		// Latitude coefficients
		double VII  = tanφ/(2*ρ*ν);  
		double VIII = tanφ/(24*ρ*ν_3)  * (5 +  3*tanφ_2 + η_2 - 9*tanφ_2*η_2);
		double IX   = tanφ/(720*ρ*ν_5) * (61 + 90*tanφ_2 + 45*tanφ_4);
		// Longitude coefficients; the OS doc uses the sec function here, but as sec = 1/cos, just use 1/cos below.
		double X    = 1.0/(cosφ*ν);
		double XI   = 1.0/(6*cosφ*ν_3) * (ν/ρ + 2*tanφ_2);
		double XII  = 1.0/(120*cosφ*ν_5) * (5 + 28*tanφ_2 + 24*tanφ_4);
		double XIIA = 1.0/(5040*cosφ*ν_7) * (61 + 662*tanφ_2 + 1320*tanφ_4 + 720*tanφ_6);

		// Apply the coefficients to powers of the difference between our easting and the projection origin easting to obtain 
		// latitude and longitude.  
		double Ediff = E - E0;
		// Latitude adustments are applied to the φ1 calculated above, the latitude value if the E value is 0. 
		double φ = φ1 - VII*pow(Ediff, 2) + VIII*pow(Ediff, 4) - IX*pow(Ediff, 6);   
		// Longitude adustments are applied to the the longitude value of the true origin.  
		double λ = λ0 + X*Ediff - XI*pow(Ediff, 3) + XII*pow(Ediff, 5) - XIIA*pow(Ediff, 7);
				
		if(m_diagnostics) {
			System.out.println();
			System.out.println("Standard values:");
			System.out.println(m_ellipsoid.asString());
			System.out.println(m_projection.asString());
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
			System.out.printf("- η²   : %e %n", η_2);
			System.out.printf("- VII  : %e %n", VII);
			System.out.printf("- VIII : %e %n", VIII);
			System.out.printf("- IX   : %e %n", IX);
			System.out.printf("- X    : %e %n", X);
			System.out.printf("- XI   : %e %n", XI);
			System.out.printf("- XII  : %e %n", XII);
			System.out.printf("- XIIA : %e %n", XIIA);
			System.out.printf("- φ    : %f rad %n", φ);
			System.out.printf("- λ    : %f rad %n", λ);
		}
		
		// Convert to our latitude and longitude objects 
		LatLong.Generator g = new LatLong.Generator(); 
		LatLong lat = g.fromRadians(LatLong.AngleType.LATITUDE, φ);
		LatLong lon = g.fromRadians(LatLong.AngleType.LONGITUDE, λ);
		LatLongPoint llp = new LatLongPoint(lat, lon);
		
		if(m_diagnostics) {
			System.out.printf("- φ    : %s %n", lat.asDMSString());
			System.out.printf("- λ    : %s %n", lon.asDMSString());			
		}
		
		return llp;
	}

	// Arc of Meridian calculation for a specified latitude. 
	// 2.21 and 2.22 in History of Retriangulation of Great Britain.
	// which references two German documents which it labels Jordan (14) and Hristow (15)
	// Perhaps also:
	// https://books.google.co.uk/books?id=f1bgBAAAQBAJ&lpg=PA300&ots=ULSkaJjjAu&dq=jordan%20eggert&pg=PA10#v=onepage&q=jordan%20eggert&f=false
	private double calculateM(double φ) {

		double n = (a-b)/(a+b);
		double n_2 = pow(n,2);
		double n_3 = pow(n,3);

		double M = b*F0 * 
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

		//Converter c = new Converter(Ellipsoid.AIRY_1830, Projection.NATIONAL_GRID);
		Converter c = new Converter(Ellipsoid.WGS_84, Projection.NATIONAL_GRID);
		c.setDiagnostics(true);
		
		if(lat.isValid() && lon.isValid()) {
			LatLongPoint pIn = new LatLongPoint(lat, lon);		
			ENPoint pOut = c.toEastingNorthing(pIn);
			System.out.println();
			System.out.println();
			System.out.println("Converted:");
			System.out.printf(" %s %n", pIn.asDMSString() );
			System.out.println("to:");
			System.out.printf(" %s %n", pOut.asString() );
		}
		else {
			System.out.println("Invalid point");
		}
		
		System.out.println();
		ENPoint p2In = new ENPoint(651409.903, 313177.270);
		
		LatLongPoint p2Out = c.toLatLong(p2In);
		System.out.println();
		System.out.println();
		System.out.println("Converted:");
		System.out.printf(" %s %n", p2In.asString() );
		System.out.println("to:");
		System.out.printf(" %s %n", p2Out.asDMSString() );
	}
}

