

pars = {
   'K': 100, # Number of breeding territories

   # Aimed for about 15 days top to end,
   'u_r':  8, # Recruitment probability sigmoid midpoint in time
   'h_r':  1/2, # Recruitment probability sigmoid slope
   'a': 3, # Maximum reproduction

   # Made this one symmetric so the optimum is at 0
   'u_s': -8, # Survival probability sigmoid midpoint in time
   'h_s': -1/2, # Survival probability sigmoid slope

   # Chose this one so its arrival date is about 3 days before the optimum (0)
   'm': 1/8, # Controls territory-competition strength: x_ss = -3.099
   's': 65/100 # Shared portion of adult and yearling survival
}
