#Written report SBL.06003 Classical models in Biology (exercises) 
#Spring semester 2025

#Student 1
#First name: Susanna
#Family name: Schärer
#Student number:

#Student 2
#First name: Elena
#Family name: Graf
#Student number:

#Student 3
#First name: Anna
#Family name: Boss
#Student number: 24-208-928

library(deSolve)
#!only de library "deSolve" is allowed!

########################
### Topic 4: Lotka-Volterra model for 2 competing species
rm(list=ls())

#R-code

lotka_volterra_competing <- function(initial_condition = c(N1 = 2.5, N2 = 0.5),
                                     parameter = c(r1 = 3, r2 = 3, alpha11 = 1, alpha12 = 0.5, alpha21 = 0.5, alpha22 = 1),
                                     scale = 0.05,
                                     equilibrium = 'eq2') {
    # Define model
    model <- function(t, y, parms) {
        with(as.list(c(y, parms)), {
            dN1 <- N1 * (r1 - alpha11 * N1 - alpha12 * N2)
            dN2 <- N2 * (r2 - alpha21 * N1 - alpha22 * N2)
            list(c(dN1, dN2))
        })
    }
    
    # Time sequence and parameters
    times <- seq(0, 30, 0.01)
    parms <- parameter
    
    # Run simulation
    out <- ode(initial_condition, times, model, parms)
    
    # Calculate equilibrium points
    # Trivial equilibria
    #eq1 <- c(0, 0)  # Both species extinct
    eq2 <- c(parms["r1"]/parms["alpha11"], 0)  # N1 only
    eq3 <- c(0, parms["r2"]/parms["alpha22"])  # N2 only
    
    # Coexistence equilibrium A*N=b -> N=A^-1*b
    A <- matrix(c(parms["alpha11"], parms["alpha12"],
                  parms["alpha21"], parms["alpha22"]), 
                nrow = 2, byrow = TRUE)
    b <- c(parms["r1"], parms["r2"])
    eq4 <- solve(A, b)  # Solve the linear system

    # Create phase plane
    plot(out[, "N1"], out[, "N2"], 
         type = "l", col = "red", lwd = 2,
         xlab = "Species 1 (N1)", ylab = "Species 2 (N2)",
         main = paste("Lotka-Volterra Competing Species\nInitial: N1 =", initial_condition["N1"], 
                      ", N2 =", initial_condition["N2"]),
         xlim = c(0, max(out[, "N1"]) * 1.1),
         ylim = c(0, max(out[, "N2"]) * 1.1))
    
    # Add isoclines
    # N1-isocline (dN1/dt = 0): N1 = 0 or r1 - alpha11*N1 - alpha12*N2 = 0
    curve((parms["r1"] - parms["alpha11"] * x)/parms["alpha12"], 
          from = 0, to = max(out[, "N1"]) * 1.1, 
          col = "blue", lwd = 2, lty = 2, add = TRUE)
    # N2-isocline (dN2/dt = 0): N2 = 0 or r2 - alpha21*N1 - alpha22*N2 = 0
    curve((parms["r2"] - parms["alpha21"] * x)/parms["alpha22"], 
          from = 0, to = max(out[, "N1"]) * 1.1, 
          col = "darkgreen", lwd = 2, lty = 2, add = TRUE)
    
    # Add equilibrium points
    #points(eq1[1], eq1[2], pch = 16, cex = 1.5, col = "black")
    if (equilibrium == "eq2"){
        points(eq2[1], eq2[2], pch = 16, cex = 1.5, col = "black")
    }
    if (equilibrium == "eq3"){
        points(eq3[1], eq3[2], pch = 16, cex = 1.5, col = "black")
    }
    if (equilibrium == "eq4"){
        points(eq4[1], eq4[2], pch = 16, cex = 1.5, col = "black")
    }
    
    
    # Add vector field
    N1_seq <- seq(0, max(out[, "N1"]) * 1.1, length.out = 15)
    N2_seq <- seq(0, max(out[, "N2"]) * 1.1, length.out = 15)
    grid <- expand.grid(N1 = N1_seq, N2 = N2_seq)
    
    for(i in 1:nrow(grid)) {
        N1 <- grid$N1[i]
        N2 <- grid$N2[i]
        derivs <- model(0, c(N1 = N1, N2 = N2), parms)[[1]]
        # Scale arrows for better visualization
        arrows(N1, N2, 
               N1 + derivs[1] * scale, 
               N2 + derivs[2] * scale, 
               length = 0.05, col = "gray")
    }
    
    # Add legend
    legend("topright", 
           legend = c("Trajectory", "N1-isocline", "N2-isocline", "Equilibrium"),
           col = c("red", "blue", "darkgreen", "black"),
           lty = c(1, 2, 2, NA), lwd = c(2, 2, 2, NA), pch = c(NA, NA, NA, 16),
           bty = "n")
    
}

# Case 1: the two species can coexist
lotka_volterra_competing(initial_condition = c(N1 = 2, N2 = 0.1), parameter = c(r1 = 1, r2 = 1, alpha11 = 0.5, alpha12 = 0.3, alpha21 = 0.2, alpha22 = 0.6), scale = 0.1, equilibrium = 'eq4')
# Case 2: species 1 outcompetes species 2
lotka_volterra_competing(initial_condition = c(N1 = 1.5, N2 = 1), parameter = c(r1 = 1, r2 = 0.8, alpha11 = 0.5, alpha12 = 0.3, alpha21 = 0.6, alpha22 = 0.4), equilibrium = 'eq3')
# Case 3: species 2 outcompetes species 1
lotka_volterra_competing(initial_condition = c(N1 = 0.2, N2 = 3), parameter = c(r1 = 3.5, r2 = 1.5, alpha11 = 0.5, alpha12 = 1, alpha21 = 1, alpha22 = 0.5), equilibrium = 'eq2')
# Case 4: priority effect: one species outcompetes the other, but it depends on the initial conditions
lotka_volterra_competing(initial_condition = c(N1 = 3.5, N2 = 3), parameter = c(r1 = 2, r2 = 2, alpha11 = 0.5, alpha12 = 1, alpha21 = 1, alpha22 = 0.5))
lotka_volterra_competing(initial_condition = c(N1 = 3, N2 = 3.5), parameter = c(r1 = 2, r2 = 2, alpha11 = 0.5, alpha12 = 1, alpha21 = 1, alpha22 = 0.5))


###Discussion of the results###
##Case 1:
#Two species can coexist when the intraspecific competition (alpha11 and alpha22) is greater than interspecific competition.
#Therefore, this following has to be true: alpha11*alpha22 > alpha12*alpha21.
#Additionally, the initial conditions has to be positive (>0) for both species.

##Case 2 and Case 3:
#For one species to win over the other species, the following two conditions have to be true (example for Species 1):
#r1/alpha12 > r2/alpha22: Species 1 can invade when Species 2 is at equilibrium and
#r2/alpha21 < r1/alpha11: Species 2 cannot invade when Species 1 is at equilibrium
#The surviving species will reach an equilibrium at carrying capacity K = r1/alpha11 or K = r2/alpha22 for Species 1 and Species 2, respectively.


#Case 4:
#


# Notes: (can be deleted)
# r:intrinsic growth rates
# alpha11, alpha22: intraspecific competition
# alpha12, alpha21: interspecific competition
# coexisting equilibrium: dN1/dt=0 and dN2/dt=0
# this leads to:
# r1 = alpha11*N1 + alpha12*N2 &
# r2 = alpha21*N1 + alpha22*N2
# matrix form: A * N = b


########################
### Topic 5: Lotka-Volterra model for 2 mutualistic species
rm(list=ls())

#R-code

lotka_volterra_mutalistic <- function(initial_condition = c(N1 = 2.5, N2 = 0.5),
                                       parameter = c(r1 = 3, r2 = 3, alpha11 = 1, gamma12 = 0.5, gamma21 = 0.5, alpha22 = 1),
                                       scale = 0.05) {
    # Define model
    model <- function(t, y, parms) {
        with(as.list(c(y, parms)), {
            dN1 <- N1 * (r1 - alpha11 * N1 + gamma12 * N2)
            dN2 <- N2 * (r2 + gamma21 * N1 - alpha22 * N2)
            list(c(dN1, dN2))
        })
    }
    
    # Time sequence and parameters
    times <- seq(0, 30, 0.01)
    parms <- parameter
    
    # Run simulation
    out <- ode(initial_condition, times, model, parms)
    
    # Calculate equilibrium points
    # Trivial equilibria
    eq1 <- c(0, 0)  # Both species extinct
    eq2 <- c(parms["r1"]/parms["alpha11"], 0)  # N1 only
    eq3 <- c(0, parms["r2"]/parms["alpha22"])  # N2 only
    
    # Coexistence equilibrium
    A <- matrix(c(parms["alpha11"], parms["alpha12"],
                  parms["alpha21"], parms["alpha22"]), 
                nrow = 2, byrow = TRUE)
    b <- c(parms["r1"], parms["r2"])
    eq4 <- solve(A, b)  # Solve the linear system
    
    # Create phase plane
    plot(out[, "N1"], out[, "N2"], 
         type = "l", col = "red", lwd = 2,
         xlab = "Species 1 (N1)", ylab = "Species 2 (N2)",
         main = paste("Lotka-Volterra Mutualistic Species\nInitial: N1 =", initial_condition["N1"], 
                      ", N2 =", initial_condition["N2"]),
         xlim = c(0, max(out[, "N1"]) * 1.1),
         ylim = c(0, max(out[, "N2"]) * 1.1))
    
    # Add isoclines
    # N1-isocline (dN1/dt = 0): N1 = 0 or r1 - alpha11*N1 - alpha12*N2 = 0
    curve((parms["r1"] - parms["alpha11"] * x)/parms["alpha12"], 
          from = 0, to = max(out[, "N1"]) * 1.1, 
          col = "blue", lwd = 2, lty = 2, add = TRUE)
    # N2-isocline (dN2/dt = 0): N2 = 0 or r2 - alpha21*N1 - alpha22*N2 = 0
    curve((parms["r2"] - parms["alpha21"] * x)/parms["alpha22"], 
          from = 0, to = max(out[, "N1"]) * 1.1, 
          col = "darkgreen", lwd = 2, lty = 2, add = TRUE)
    
    # Add equilibrium points
    #points(eq1[1], eq1[2], pch = 16, cex = 1.5, col = "black")
    points(eq2[1], eq2[2], pch = 16, cex = 1.5, col = "black")
    points(eq3[1], eq3[2], pch = 16, cex = 1.5, col = "black")
    points(eq4[1], eq4[2], pch = 16, cex = 1.5, col = "black")
    
    # Add vector field
    N1_seq <- seq(0, max(out[, "N1"]) * 1.1, length.out = 15)
    N2_seq <- seq(0, max(out[, "N2"]) * 1.1, length.out = 15)
    grid <- expand.grid(N1 = N1_seq, N2 = N2_seq)
    
    for(i in 1:nrow(grid)) {
        N1 <- grid$N1[i]
        N2 <- grid$N2[i]
        derivs <- model(0, c(N1 = N1, N2 = N2), parms)[[1]]
        # Scale arrows for better visualization
        arrows(N1, N2, 
               N1 + derivs[1] * scale, 
               N2 + derivs[2] * scale, 
               length = 0.05, col = "gray")
    }
    
    # Add legend
    legend("topright", 
           legend = c("Trajectory", "N1-isocline", "N2-isocline", "Equilibrium"),
           col = c("red", "blue", "darkgreen", "black"),
           lty = c(1, 2, 2, NA), lwd = c(2, 2, 2, NA), pch = c(NA, NA, NA, 16),
           bty = "n")
    
}



# (1) there exists a positive equilibrium points which is stable
lotka_volterra_mutalistic(initial_condition =  c(N1 = 0.5, N2 = 0.5), parameter = c(r1=1, r2=0.5, alpha11=2, gamma12=1.5, gamma21=1.11, alpha22=1.5))
# (2)  there exists a positive equilibrium points, but which is unstable
lotka_volterra_mutalistic(initial_condition =  c(N1 = 0.9, N2 = 0.9), parameter = c(r1=0.5, r2=-1, alpha11=1.2, gamma12=0.7, gamma21=0.7, alpha22=0.2))

#Discussion of the results

