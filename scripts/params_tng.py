basepath    = '/home/cemoser/TNG100-3/output'
serial      = True
search_radius = 50 #10 scaled
lims        = [10e-2,10e4]  #[0.01,10] # scaled radius
bins        = 25
mass_option = 1 #Option 1 = upper/lower mass bounds, option 2 = central mass plus percentage range
h=0.677
mass_low    = h*10**9.5 #10**9.5 #10.**11 #these are the halo limits 
mass_high   = h*10**12.3 #10**12.3 #5.*10**14
mass_option_percent = 0.2 #This is +- percentage of mass option 2, as a fraction 
volweight   = True #To convert mass to density
scaled_radius = False 
mass_kind   = 'stellar' #search by halo or stellar mass, options='stellar','halo' 
