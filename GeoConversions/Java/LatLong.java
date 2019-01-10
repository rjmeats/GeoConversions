
public class LatLong {

	// Enum used for distinguishing latitude and longitudes
	
	public static enum AngleType {
		
		LATITUDE(90, "N", "S"), LONGITUDE(180, "E", "W");
		
		private double m_maxDegrees;
		private String m_positiveDirectionIndicator;
		private String m_negativeDirectionIndicator;
		
		AngleType(double maxDegrees, String pos, String neg) {
			m_maxDegrees = maxDegrees;
			m_positiveDirectionIndicator = pos;
			m_negativeDirectionIndicator = neg;
		}
		
		public boolean inDegreesRange(double degrees) {
			return Math.abs(degrees) <= m_maxDegrees;
		}

		public String getPositiveDirectionString() {
			return m_positiveDirectionIndicator;
		}
		
		public String getNegativeDirectionString() {
			return m_negativeDirectionIndicator;
		}		

	}

	// Main latitude/longitude value class
	AngleType m_type;
	PositionAngle m_angle;
	String m_direction;
	private boolean m_valid;
	
	protected LatLong(AngleType type) {
		m_type = type;
		m_angle = PositionAngle.fromDegrees(0);
		m_direction = "";
		m_valid = true;
	}

	
	public String asString() {
		String o = PositionAngle.DMS_Angle.DEGREE_SYMBOL;		
		if(isValid()) {
			return String.format("%.4f%s %s", Math.abs(m_angle.degrees()), o, (m_direction.length() > 0 ? m_direction : "0"));
		}
		else {
			return "*** Not valid ***";
		}
	}

	public String asDMSString() {
		if(isValid()) {
			return String.format("%s %s", m_angle.asUnsignedDMSString(), (m_direction.length() > 0 ? m_direction : "0"));
		}
		else {
			return "*** Not valid ***";
		}		
	}

	public String asDetailString() {
		if(isValid()) {
			return String.format("%s : %s", m_angle.asString(), (m_direction.length() > 0 ? m_direction : "0"));
		}
		else {
			return "*** Not valid ***";
		}
	}
	
	public boolean isValid() { return m_valid; }
	public double asDegrees() { return m_angle.m_degrees; }
	public double asRadians() { return m_angle.m_radians; }

	// =====================================================================================
	
	// Class for generating latitude / longitude objects
	
	public static class Generator {
		
		public LatLong fromDegrees(AngleType type, double degrees) {
			LatLong l = new LatLong(type);
			PositionAngle a = PositionAngle.fromDegrees(degrees);
			String direction = "";
			if(a.isValid()) {
				if(degrees > 0) {
					direction = type.getPositiveDirectionString();
				}
				else if(degrees < 0) {
					direction = type.getNegativeDirectionString();
				}
				else if(degrees == 0.0) {
					direction = "";
				}
				
				l.m_valid = true;
				l.m_angle = a;
				l.m_direction = direction;
			}
			else {
				l.m_valid = false;
			}
	
			if(l.isValid()) {
				if(!type.inDegreesRange(l.m_angle.degrees())) {
					l.m_valid = false;
				}
			}
			
			return l;
		}
	
		public LatLong fromRadians(AngleType type, double radians) {
			double degrees = Math.toDegrees(radians);
			return fromDegrees(type, degrees);
		}
	
		public LatLong fromDMS(AngleType type, int degrees, int minutes, double seconds, String direction) {
			LatLong l = new LatLong(type);
			int sign = 0;
			boolean directionOK = true;
			String myDirection = "";
			if(direction.toUpperCase().equals(type.getPositiveDirectionString()) || direction.equals("+")) {
				sign = +1;
				myDirection = type.getPositiveDirectionString();
			}
			else if(direction.toUpperCase().equals(type.getNegativeDirectionString()) || direction.equals("-")) {
				sign = -1;
				myDirection = type.getNegativeDirectionString();
			}
			else if(direction.toUpperCase().equals("")) {
				sign = 0;
				myDirection = "";
			}
			else {
				directionOK = false;
			}
			
			if(directionOK) {
				PositionAngle a = PositionAngle.fromDMS(degrees, minutes, seconds, sign);
				if(a.isValid()) {
					l.m_valid = true;
					l.m_angle = a;
					l.m_direction = myDirection; 
				}
				else {
					l.m_valid = false;
				}
			}
			else {
				l.m_valid = false;
			}				
			
			if(l.isValid()) {
				if(!type.inDegreesRange(l.m_angle.degrees())) {
					l.m_valid = false;
				}
			}
			
			return l;
		}
		
	}	
	
	public static void main(String args[]) {
		
		LatLong.Generator g = new LatLong.Generator(); 
		LatLong lat = g.fromDegrees(AngleType.LATITUDE, 52.0);
		System.out.println(lat.asString());
		
		lat = g.fromDegrees(AngleType.LATITUDE, -10.2);
		System.out.println(lat.asString());
		
		lat = g.fromDMS(AngleType.LATITUDE, 52, 10, 30, "s");
		System.out.println(lat.asString());
	}
}
