###Calculates theoretical expecations for splashback radii of the clusters
from colossus.cosmology import cosmology
from colossus.halo import concentration
from colossus.halo.mass_defs import changeMassDefinition
import colossus.halo.splashback as sb

h = 0.7 #Hubble parameter
cosmology.setCosmology('bolshoi')


def calc_stat(mass, redshift, mass_definition = "200m"):
    if mass_definition == "200m":
        M200m = mass
        c200m = concentration.concentration(M200m, mass_definition, redshift)
        M500c, R500c, c500c = changeMassDefinition(M200m, c200m, redshift,'200m', '500c', profile = 'nfw')
    elif mass_definition == "500c":
        M500c = mass
	c500c = concentration.concentration(M500c, '500c', redshift)
        M200m, R200m, c200m = changeMassDefinition(M500c, c500c, redshift,'500c', '200m')

    Rsp = sb.splashbackRadius(redshift, '200m', M = M200m)[0]*(1 + redshift)
    return M200m, M500c, R200m, Rsp


if __name__=="__main__":
   
    #RedMaPPer
    mass = 1.8e14
    redshift = 0.24
    M200m, M500c, R200m = calc_stat(mass, redshift, mass_defintion)
    print("RedMaPPer:")
    print("######################################")
    print("R200m is: %E kpc/h" % R200m)
    print("M500c is %E Msol/h and M200m is %E Msol/h" % (M500c, M200m))
    print("The splashback radius is: %f kpc/h" % Rsp)
    print("######################################")

    #SZ
    mass = 4.3e14*h
    redshift = 0.177
    M200m, M500c, R200m = calc_stat(mass, redshift, mass_defintion)
    print("SZ:")
    print("######################################")
    print("R200m is: %E kpc/h" % R200m)
    print("M500c is %E Msol/h and M200m is %E Msol/h" % (M500c, M200m))
    print("The splashback radius is: %f kpc/h" % Rsp)
    print("######################################")

    #XR
    mass = 2.3e14*h
    redshift = 0.186
    M200m, M500c, R200m = calc_stat(mass, redshift, mass_defintion)
    print("XR:")
    print("######################################")
    print("R200m is: %E kpc/h" % R200m)
    print("M500c is %E Msol/h and M200m is %E Msol/h" % (M500c, M200m))
    print("The splashback radius is: %f kpc/h" % Rsp)
    print("######################################")
