/*
 * Convert latitude and longitude to OS National Grid Easting and Northing.
 * 
 * Using formulas from Annexe C of 'A Guide to the Coordinate Systems in Great Britain',
 * which is available at:
 * 
 * 	https://www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf
 * 
 * NB Longitude east of meridian is postive.
 */

import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.tan;

public class Converter {

	// Airy 1830 ellipsoid values, from Annexe A 
	public static double SEMI_MAJOR_AXIS = 6_377_563.396;
	public static double SEMI_MINOR_AXIS = 6_356_256.909;

	public static double EASTING_OF_TRUE_ORIGIN = 400_000;
	public static double NORTHING_OF_TRUE_ORIGIN = -100_000;
	public static double SCALE_FACTOR_ON_CENTRAL_MERIDIAN = 0.9996012717;
	public static double LATITUDE_OF_TRUE_ORIGIN_DEGREES = 49.0;	// Degrees N
	public static double LONGITUDE_OF_TRUE_ORIGIN_DEGREES = -2.0;   // Degrees W
	public static double LATITUDE_OF_TRUE_ORIGIN = degreesToRadians(LATITUDE_OF_TRUE_ORIGIN_DEGREES, 0, 0); 
	public static double LONGITUDE_OF_TRUE_ORIGIN = degreesToRadians(LONGITUDE_OF_TRUE_ORIGIN_DEGREES, 0, 0); 
	
	public static double degreesToRadians(double degrees, double minutes, double seconds) {
		double decimalDegrees = degrees + minutes/60.0 + seconds/3600.0;
		return Math.toRadians(decimalDegrees);
	}
	
	public static void toEastingsNorthings(double latitude, double longitude) {
		
		// Use same symbols as in Annexe C 
		double a = SEMI_MAJOR_AXIS;
		double b = SEMI_MINOR_AXIS;
		
		// From chapter B.1
		double e_2 = ((a*a) - (b*b)) / (a*a);
		
		// From chapter A.2
		double E0 = EASTING_OF_TRUE_ORIGIN;
		double N0 = NORTHING_OF_TRUE_ORIGIN;
		double F0 = SCALE_FACTOR_ON_CENTRAL_MERIDIAN;
		
		double φ0 = LATITUDE_OF_TRUE_ORIGIN;
		double λ0 = LONGITUDE_OF_TRUE_ORIGIN;

		double φ = latitude;
		double λ = longitude;

		double n = (a-b)/(a+b);
		double n_2 = n*n;
		double n_3 = n*n*n;
		
		double sinφ = sin(φ);
		double sinφ_2 = pow(sin(φ), 2);
		double cosφ = cos(φ);
		double cosφ_3 = pow(cos(φ), 3);
		double cosφ_5 = pow(cos(φ), 5);
		double tanφ_2 = pow(tan(φ), 2);
		double tanφ_4 = pow(tan(φ), 4);
		
		double ν = a * F0 * Math.pow(1 - (e_2 * sinφ_2), -0.5); 
		double ρ = a * F0 * (1 - e_2) * Math.pow(1 - (e_2 * sinφ_2), -1.5);		
		double ηη = (ν / ρ) - 1;  
		
		double M = b * F0 * 
					(   (  ( 1 + n + (5/4)*n_2 +  (5/4)*n_3 ) * (φ-φ0)  )                             -
					    (  (   3*n +     3*n_2 + (21/8)*n_3 ) * sin(φ-φ0)     * cos(φ+φ0) )           + 
					    (  (        (15/8)*n_2 + (15/8)*n_3 ) * sin(2*(φ-φ0)) * cos(2*(φ+φ0)) )       -
					    (  (                    (35/24)*n_3 ) * sin(3*(φ-φ0)) * cos(3*(φ+φ0)) )
					);

		double I    = M + N0;
		double II   = (ν/2)   * sinφ * cosφ;
		double III  = (ν/24)  * sinφ * cosφ_3 * (5 -     tanφ_2 + 9*ηη );
		double IIIA = (ν/720) * sinφ * cosφ_5 * (61 - 58*tanφ_2 + tanφ_4);
		double IV   = ν       * cosφ;
		double V    = (ν/6)   * cosφ_3 * (ν/ρ -  tanφ_2);
		double VI   = (ν/120) * cosφ_5 * (5 - 18*tanφ_2 + tanφ_4 + 14*ηη - 58*tanφ_2*ηη); 

		double λDiff = λ-λ0;
		double N = I  + II * pow(λDiff, 2) + III * pow(λDiff, 4) + IIIA * pow(λDiff, 6);
		double E = E0 + IV * pow(λDiff, 1) + V   * pow(λDiff, 3) + VI   * pow(λDiff, 5);
		
		System.out.println("Standard values:");
		System.out.printf("- a    : %f %n", a);
		System.out.printf("- b    : %f %n", b);
		System.out.printf("- e*e  : %f %n", e_2);
		System.out.println();
		System.out.printf("- E0   : %f %n", E0);
		System.out.printf("- N0   : %f %n", N0);
		System.out.printf("- F0   : %f %n", F0);
		System.out.println();
		System.out.printf("- φ0   : %f %n", φ0);
		System.out.printf("- λ0   : %f %n", λ0);
		System.out.println();
		System.out.printf("- n    : %f %n", n);
		System.out.println();
		System.out.println("Values for the specific latitude/longitude:");
		System.out.println();
		System.out.printf("- φ    : %f %n", φ);
		System.out.printf("- λ    : %f %n", λ);
		System.out.println();
		System.out.printf("- ν    : %f %n", ν);
		System.out.printf("- ρ    : %f %n", ρ);
		System.out.printf("- ηη   : %f %n", ηη);
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
	}
	
	public static void main(String args[]) {
				
		toEastingsNorthings(degreesToRadians(52, 39, 27.2531), degreesToRadians(1, 43, 4.5177));
	}
}

