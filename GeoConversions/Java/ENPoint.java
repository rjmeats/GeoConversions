
public class ENPoint {

	// Full values
	double m_easting;
	double m_northing;
	
	ENPoint(double e, double n) {
		m_easting = e;
		m_northing = n;
	}
	
	public String asString() {
		return String.format("E=%.3f / N=%.3f", m_easting, m_northing);
	}
}
