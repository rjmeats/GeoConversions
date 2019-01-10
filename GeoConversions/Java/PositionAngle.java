import java.lang.Math;

public class PositionAngle {

	public boolean m_valid;
	public double m_degrees;		
	public double m_radians;	
	public DMS_Angle m_dms;	

	public static PositionAngle fromDegrees(double degrees) {
		PositionAngle a = new PositionAngle();		
		if(Math.abs(degrees) < 360.0) {
			a.m_valid = true;
			a.m_degrees = degrees;
			a.m_radians = Math.toRadians(a.m_degrees);
			a.m_dms = DMS_Angle.fromDegrees(a.m_degrees); 
		}
		else {
			a.m_valid = false;
		}
		return a;
	}
	
	public static PositionAngle fromRadians(double radians) {
		PositionAngle a = new PositionAngle();		
		if(Math.abs(radians) < 2 * Math.PI) {
			a.m_valid = true;
			a.m_radians = radians;
			a.m_degrees = Math.toDegrees(a.m_degrees);
			a.m_dms = DMS_Angle.fromDegrees(a.m_degrees); 
		}
		else {
			a.m_valid = false;
		}
		return a;
	}
	
	public static PositionAngle fromDMS(int degrees, int minutes, double seconds, int sign) {
		PositionAngle a = new PositionAngle();		
		DMS_Angle dms = DMS_Angle.fromDMS(degrees, minutes, seconds, sign);
		if(dms.m_valid) {
			a = fromDMS(dms);
		}
		else {
			a = new PositionAngle();
			a.m_valid = false;
		}
		return a;
	}
	
	public static PositionAngle fromDMS(DMS_Angle dms) {
		PositionAngle a = new PositionAngle();
		
		if(dms.m_valid) {
			a.m_valid = true;
			a.m_dms = DMS_Angle.fromDMS(dms);
			double unsignedDecimalDegrees = dms.m_degrees + dms.m_minutes/60.0 + dms.m_seconds/3600.0;
			a.m_degrees = dms.m_sign * unsignedDecimalDegrees;
			a.m_radians = Math.toRadians(a.m_degrees);
		}
		else {
			a.m_valid = false;
		}
		
		return a;
	}
	
	private PositionAngle() {
		m_valid = true;
		m_degrees = 0.0;		
		m_radians = 0.0;	
		m_dms = new DMS_Angle();
	}

	public boolean isValid() { return m_valid; }
	public double degrees() { return m_degrees; }
	public double radians() { return m_radians; }
	public DMS_Angle dms() { return m_dms; }
	
	public String asString() {
		return isValid() ? 
				String.format("%f%s = %f radian%s = %s", 
								degrees(), DMS_Angle.DEGREE_SYMBOL, 
								radians(), radians() == 1 ? "" : "s",
								m_dms.asString())
				: "*** Invalid angle ***";
	}

	public String asDMSString() {
		return m_dms.asString();
	}

	public String asUnsignedDMSString() {
		return m_dms.asUnsignedString();
	}

	// ==============================================================================================
	
	// Represent angle as degrees, minutes and seconds, with a sign (values are stored as positive)
	// NB Some issues with roundings
	public static class DMS_Angle { 
		
		public static final String DEGREE_SYMBOL = "°";
		public static final String MINUTE_SYMBOL = "′";
		public static final String SECOND_SYMBOL = "″";
		
		private boolean m_valid;
		private int m_sign;			// +1 or 0 or -1 
		private int m_degrees;		// 0 - 359
		private int m_minutes;		// 0 - 59
		private double m_seconds;	// 0.0 - 59.99....

		public static DMS_Angle fromDMS(int degrees, int minutes, double seconds, int sign) {
			DMS_Angle a = new DMS_Angle();

			// Check for valid combinations
			if(
					(sign == 0 && degrees == 0 && minutes == 0 && seconds == 0.0) ||
					(Math.abs(sign) == 1 && 
						degrees >= 0 && degrees < 360 && 
						minutes >= 0 && minutes < 60 &&
						seconds >= 0.0 && seconds <= 60.0)		// 59.99999999999999 gets rounded to 60 so allow ????
					)
			{
				a.m_valid = true;
				a.m_sign = sign;
				a.m_degrees = degrees;
				a.m_minutes = minutes;
				a.m_seconds = seconds;
				
				// If zero, make sure sign is 0
				if(a.m_degrees == 0 && a.m_minutes == 0 && a.m_seconds == 0.0) {
					a.m_sign = 0;
				}
			}
			else {
				a.m_valid = false;
			}
			return a;
		}
		
		public static DMS_Angle fromDegrees(double degrees) {
			DMS_Angle a = new DMS_Angle();			
			if(Math.abs(degrees) < 360.0) {
				a.m_sign = PositionAngle.getSign(degrees);
				if(a.m_sign == 0) {
					// Ensure exactly zero.
					a.m_degrees = 0;
					a.m_minutes = 0;
					a.m_seconds = 0.0;
				}
				else {
					// Work out degrees and minutes as integers, retaining decimal places for the seconds
					double absDegrees = Math.abs(degrees);
					a.m_degrees = (int) absDegrees;
					double fractionOfADegree = absDegrees - a.m_degrees*1.0;
					double minutes = fractionOfADegree*60.0;
					a.m_minutes = Math.abs((int) (minutes+0.00000001));		// Make e.g. 10.2 as 10d 12m 0s instead of 10d 11m 60s ????
					double fractionOfAMinute = minutes - a.m_minutes;
					double seconds = Math.abs(fractionOfAMinute*60.0);
					a.m_seconds = seconds;
				}
			}
			else {
				a.m_valid = false;
			}
			return a;
		}
		
		// Copy an existing object
		public static DMS_Angle fromDMS(DMS_Angle dms) {
			DMS_Angle a = new DMS_Angle();
			a.m_valid = dms.m_valid;
			a.m_sign = dms.m_sign;
			a.m_degrees = dms.m_degrees;
			a.m_minutes = dms.m_minutes;
			a.m_seconds = dms.m_seconds;
			return a;
		}

		private DMS_Angle() {
			m_valid = true;
			m_degrees = 0;
			m_minutes = 0;
			m_seconds = 0.0;
			m_sign = 0;
		}
		
		public boolean isValid() { return m_valid; }
		public int sign() { return m_sign; }
		public int d() { return m_degrees; }
		public int m() { return m_minutes; }
		public double s() { return m_seconds; }
		
		public String asString() {
			String signString = "?";
			switch(sign()) {
				case -1 : signString = "-"; break;
				case  0 : signString = "0";  break;
				case +1 : signString = "+"; break;
			}
			
			return isValid() ? 
					String.format("%d%s %02d%s %f%s (%s)", 
									d(), DEGREE_SYMBOL, 
									m(), MINUTE_SYMBOL,
									s(), SECOND_SYMBOL,
									signString) 
					: "*** Invalid DMS ***";
		}
		
		public String asUnsignedString() {
			return isValid() ? 
					String.format("%d%s %02d%s %f%s", 
									d(), DEGREE_SYMBOL, 
									m(), MINUTE_SYMBOL,
									s(), SECOND_SYMBOL)
					: "*** Invalid DMS ***";
		}
	}
	
	// The Java Math library doesn't seem to have a method to return an int 1/0/-1 value according to the sign of a number, so 
	// we need to implement them here.
	
	private static int getSign(double n) {
		if(n < 0.0) return -1;
		else if(n == 0.0)  return 0;
		else return +1;
	}
	
	// Some tests
	public static void main(String args[]) {
		
		DMS_Angle dms = DMS_Angle.fromDegrees(10.2);
		System.out.println(dms.asString());

		PositionAngle pa = PositionAngle.fromDegrees(360/Math.PI);
		System.out.println(pa.asString());
		
		pa = PositionAngle.fromDMS(10, 11, 59.99999990, 1);
		System.out.println(pa.asString());
	}
}
