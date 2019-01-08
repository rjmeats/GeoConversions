
public class LatLongPoint {

	public LatLong m_latitude;
	public LatLong m_longitude;
	
	public LatLongPoint(LatLong lat, LatLong lon) {
		m_latitude = lat;
		m_longitude = lon;
	}
	
	public String asString() {
		return String.format("Lat=%s / Long=%s", m_latitude.asString(), m_longitude.asString());
	}
}
