
public class ConverterLambda {

	public HandlerResponse handler(HandlerInput input) {

		Converter con = new Converter(Converter.Ellipsoid.AIRY_1830, Converter.Projection.NATIONAL_GRID);
		con.setDiagnostics(true);

		double dLat = Double.parseDouble(input.getLatitude());
		double dLon = Double.parseDouble(input.getLongitude());

		LatLong.Generator g = new LatLong.Generator(); 
		LatLong lat = g.fromDegrees(LatLong.AngleType.LATITUDE, dLat);
		LatLong lon = g.fromDegrees(LatLong.AngleType.LONGITUDE, dLon);

		ENPoint p = new ENPoint(0, 0);
		if(lat.isValid() && lon.isValid()) {
			LatLongPoint pIn = new LatLongPoint(lat, lon);		
			ENPoint pOut = con.toEastingNorthing(pIn);
			System.out.println();
			System.out.println();
			System.out.println("Converted:");
			System.out.printf(" %s %n", pIn.asDMSString() );
			System.out.println("to:");
			System.out.printf(" %s %n", pOut.asString() );
			p = pOut;			
		}
		else {
			System.out.println("Invalid point");
		}
		
		return new HandlerResponse(p);
	}
	
	// Input is JSON providing latitude and longitude
	public static class HandlerInput {
		private String m_latitude;
		private String m_longitude;
		
		public String getLatitude() { return m_latitude; }
		public void setLatitude(String latitude) { m_latitude = latitude; }
		
		public String getLongitude() { return m_longitude; }
		public void setLongitude(String longitude) { m_longitude = longitude; } 
	}
	
	public static class HandlerResponse {
		private ENPoint m_point;
		HandlerResponse(ENPoint point) { this.m_point = point; }
		
		public int getE() { return (int)Math.round(m_point.m_easting); }
		public int getN() { return (int)Math.round(m_point.m_northing); }
	}
}
