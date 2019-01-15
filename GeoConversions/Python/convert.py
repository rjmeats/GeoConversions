# Convert between Latitude/Longitude and National Grid coordindates

import numpy as np
import math

def toDMS(degrees) :

    sign = int(np.sign(degrees))
    abs_degrees = abs(degrees)
    i_degrees = int(abs_degrees)
    i_minutes = int((abs_degrees - i_degrees) * 60)
    seconds =   (abs_degrees - i_degrees - i_minutes/60.0) * 3600

    return (sign, i_degrees, i_minutes, seconds)

def toDegreesFromDMS(sign, i_degrees, i_minutes, seconds) :

    return sign * (i_degrees + i_minutes/60 + seconds/60/60)

def toRadians(degrees) :
    return math.radians(degrees)

def toDegreesFromRadians(radians) :
    return math.degrees(radians)

def convertToEastingNorthing(ellipsoid, projection, lat_d, long_d) :

    print()
    print("Values derived for the specific latitude/longitude value:")

    a = ellipsoid[0]
    b = ellipsoid[1]
    e_2 = ((a*a) - (b*b)) / (a*a)

    F0 = projection[0]
    #φ0 = projection[1]
    λ0 = projection[2]
    E0 = projection[3]
    N0 = projection[4]

    φ = toRadians(lat_d)                    # phi = latitude in radians
    λ = toRadians(long_d)		            # lambda = longitude in radians 

    # Pre-calculate some trig functions and their powers of the latitude value
    sinφ = math.sin(φ)
    cosφ = math.cos(φ)
    tanφ = math.tan(φ)
    sinφ_2 = pow(sinφ, 2)
    cosφ_3 = pow(cosφ, 3)
    cosφ_5 = pow(cosφ, 5)
    tanφ_2 = pow(tanφ, 2)
    tanφ_4 = pow(tanφ, 4)
		
    # Some intermediate values, using variable names matching OS document. Not clear what they represent! 
    # M - something relating to the Meridian arc (northing). Complex calculation also used in converting in the other direction, so has its own method
    # ν = nu
    # ρ = rho
    # η = eta
    M = calculateM(ellipsoid, projection, φ)
    ν = a*F0 * pow(1 - (e_2 * sinφ_2), -0.5)
    ρ = a*F0 * (1 - e_2) * pow(1 - (e_2 * sinφ_2), -1.5)
    η_2 = (ν/ρ) - 1
	
    # Coefficients for components of the final N and E calculation. Again using OS naming.
    # Northing coefficients
    I    = M		# OS doc adds N0 here, but I do it below to make N and E calculations more similar.
    II   = (ν/2)   * sinφ * cosφ
    III  = (ν/24)  * sinφ * cosφ_3 * (5 -     tanφ_2 + 9*η_2 )
    IIIA = (ν/720) * sinφ * cosφ_5 * (61 - 58*tanφ_2 + tanφ_4)
    # Easting coefficients
    IV   = ν       * cosφ
    V    = (ν/6)   * cosφ_3 * (ν/ρ -  tanφ_2)
    VI   = (ν/120) * cosφ_5 * (5 - 18*tanφ_2 + tanφ_4 + 14*η_2 - 58*tanφ_2*η_2)

    # Apply the coefficients to powers of the difference between our longitude and the projection origin longitude to obtain N and E values  
    λDiff = λ-λ0
    N = N0 + I + II*pow(λDiff, 2) + III*pow(λDiff, 4) + IIIA*pow(λDiff, 6)
    E = E0 +     IV*pow(λDiff, 1) + V*pow(λDiff, 3)   + VI*pow(λDiff, 5)

    print()
    print("- φ    : ", φ)
    print("- λ    : ", λ)
    print()
    print("- ν    : ", ν)
    print("- ρ    : ", ρ)
    print("- η²   : ", η_2)
    print("- M    : ", M)
    print("- I    : ", I+N0)		# Add N0 here to match OS calculation value
    print("- II   : ", II)
    print("- III  : ", III)
    print("- IIIA : ", IIIA)
    print("- IV   : ", IV)
    print("- V    : ", V)
    print("- VI   : ", VI)
    print()
    print("- N   : ", N)
    print("- E   : ", E)

    return (E, N)

def convertToLatLong(ellipsoid, projection, E, N) :

    print()
    print("Values derived for the specific N/E coordinates:")
    
    a = ellipsoid[0]
    b = ellipsoid[1]
    e_2 = ((a*a) - (b*b)) / (a*a)

    F0 = projection[0]
    φ0 = projection[1]
    λ0 = projection[2]
    E0 = projection[3]
    N0 = projection[4]

    # No exact formula, use an iterative approach to calculate what the latitude would be (φ1) for this
    # northing if the easting was that of the true origin (i.e. E = E0). Calculate the 'M' value for
    # the latitude until M is very close to the N-N0 northing offset.

    diffLimit = 0.001 * 0.01	# Keep iterating until difference is 0.01 mm or less		
    M = 0
    φ1 = φ0						# Start with the latitude set to that of the true origin. 
    loops = 0
    print()
    while True :
        loops += 1
        old_φ = φ1
        φ1 = (N-N0-M)/(a*F0) + old_φ
        M = calculateM(ellipsoid, projection, φ1)
        
        diffCheck = abs(N-N0-M)
        print("Loop ", loops, ": new_φ = ", φ1, "M=", M, "diff=", diffCheck)
        
        if diffCheck < diffLimit :
            break
        
        # Protect against non-convergence
        if loops > 100 :
            print("Excessive loops calculating latitude")
            break

    # Pre-calculate some trig functions and their powers of the latitude value
    sinφ = math.sin(φ1)
    cosφ = math.cos(φ1)
    tanφ = math.tan(φ1)
    sinφ_2 = pow(sinφ,2)
    tanφ_2 = pow(tanφ,2)
    tanφ_4 = pow(tanφ,4)
    tanφ_6 = pow(tanφ,6)

    # As before, some intermediate values, using variable names matching OS document. Not clear what they represent! 
    # ν = nu
    # ρ = rho
    # η = eta
    ν = a*F0 * pow(1 - (e_2 * sinφ_2), -0.5) 
    ρ = a*F0 * (1 - e_2) * pow(1 - (e_2 * sinφ_2), -1.5)		
    η_2 = (ν/ρ) - 1  
    ν_3 = pow(ν,3)
    ν_5 = pow(ν,5)
    ν_7 = pow(ν,7)
    
    # Coefficients for components of the final lat and long calculation. Again using OS naming.
    # Latitude coefficients
    VII  = tanφ/(2*ρ*ν)
    VIII = tanφ/(24*ρ*ν_3)  * (5 +  3*tanφ_2 + η_2 - 9*tanφ_2*η_2)
    IX   = tanφ/(720*ρ*ν_5) * (61 + 90*tanφ_2 + 45*tanφ_4)
    # Longitude coefficients the OS doc uses the sec function here, but as sec = 1/cos, just use 1/cos below.
    X    = 1.0/(cosφ*ν)
    XI   = 1.0/(6*cosφ*ν_3) * (ν/ρ + 2*tanφ_2)
    XII  = 1.0/(120*cosφ*ν_5) * (5 + 28*tanφ_2 + 24*tanφ_4)
    XIIA = 1.0/(5040*cosφ*ν_7) * (61 + 662*tanφ_2 + 1320*tanφ_4 + 720*tanφ_6)

    # Apply the coefficients to powers of the difference between our easting and the projection origin easting to obtain 
    # latitude and longitude.  
    Ediff = E - E0
    # Latitude adustments are applied to the φ1 calculated above, the latitude value if the E value is 0. 
    φ = φ1 - VII*pow(Ediff, 2) + VIII*pow(Ediff, 4) - IX*pow(Ediff, 6)
    # Longitude adustments are applied to the the longitude value of the true origin.  
    λ = λ0 + X*Ediff - XI*pow(Ediff, 3) + XII*pow(Ediff, 5) - XIIA*pow(Ediff, 7)

    φ_deg = toDegreesFromRadians(φ)
    λ_deg = toDegreesFromRadians(λ)

    print()
    print("- E    : ", E)
    print("- N    : ", N)
    print()
    print("- φ1   : ", φ1)
    print()
    print("- ν    : ", ν)
    print("- ρ    : ", ρ)
    print("- η²   : ", η_2)
    print("- VII  : ", VII)
    print("- VIII : ", VIII)
    print("- IX   : ", IX)
    print("- X    : ", X)
    print("- XI   : ", XI)
    print("- XII  : ", XII)
    print("- XIIA: ", XIIA)
    print()
    print("- φ rad : ", φ)
    print("- λ rad : ", λ)
    print("- φ deg : ", φ_deg)
    print("- λ deg : ", λ_deg)
    print("- φ     : ", toDMS(φ_deg))
    print("- λ     : ", toDMS(λ_deg))

    return (φ_deg, λ_deg)

def calculateM(ellipsoid, projection, φ) :

    a = ellipsoid[0]
    b = ellipsoid[1]
    n = (a-b)/(a+b)
    n_2 = pow(n,2)
    n_3 = pow(n,3)

    F0 = projection[0]
    φ0 = projection[1]

    M = b*F0 *                                                                                  \
        (   (  ( 1 + n + (5/4)*n_2 +  (5/4)*n_3 ) * (φ-φ0)  )                             -     \
            (  (   3*n +     3*n_2 + (21/8)*n_3 ) * math.sin(φ-φ0)     * math.cos(φ+φ0) )           +     \
            (  (        (15/8)*n_2 + (15/8)*n_3 ) * math.sin(2*(φ-φ0)) * math.cos(2*(φ+φ0)) )       -     \
            (  (                    (35/24)*n_3 ) * math.sin(3*(φ-φ0)) * math.cos(3*(φ+φ0)) )             \
        )

    return M

airy1830Ellipsoid = (6377563.396, 6356256.909)
nationalGridProjection = (0.9996012717, toRadians(toDegreesFromDMS(+1, 49, 0, 0)), toRadians(toDegreesFromDMS(-1, 2, 0, 0)), 400000, -100000)

if __name__ == "__main__" :

    lat_d = toDegreesFromDMS(1, 52, 39, 27.2531)
    lon_d = toDegreesFromDMS(1, 1, 43, 4.5177)
    ellipsoid = airy1830Ellipsoid
    projection = nationalGridProjection

    convertToEastingNorthing(ellipsoid, nationalGridProjection, lat_d, lon_d)

    print()
    print("=========================================================")
    print()

    e = 651409.903
    n = 313177.270
    convertToLatLong(ellipsoid, projection, e, n)
