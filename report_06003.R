#Written report SBL.06003 Classical models in Biology (exercises)
#Spring semester 2025

#Student 1
#First name: Susanna
#Family name: Schärer
#Student number: 19-811-652

#Student 2
#First name: Elena
#Family name: Graf
#Student number:20-701-132

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


#############################################################


### Topic 5: Lotka-Volterra model for 2 mutualistic species

# We describe the interaction between two species (N1, N2) given by the following equations :
# dN1/dt = N1 * (r1 - alpha11*N1 + gamma12*N2)
# dN2/dt = N2 * (r2 + gamma21*N1 - alpha22*N2)
# Parameters
# r1, r_2 : intrinsic growth rates (can be positive or negative)
# alpha11, alpha22: intra-specific competition (negative self-interaction), both > 0
# gamma12, gamma21 : inter-specific competition (positive inter-specific interaction), both > 0

rm(list=ls())

#R-code

lotka_volterra_mutalistic <- function(initial_condition = c(N1 = 2.5, N2 = 0.5),
                                      parameters = c(r1 = 3, r2 = 3, alpha11 = 1, gamma12 = 0.5, gamma21 = 0.5, alpha22 = 1),
                                      scale = 0.005) {
  # initial_condition : starting population sizes for species 1 and 2
  # parameters : vector containing the six model parameters (alpha11/22, gamma12/21, r1/2)
  # scale: a value that determines the length of arrows in the vector field

  # ------------- Simulation ---------------------------------------------------

  # We define the model using the provided coupled ODEs
  # The function returns the instantaneous rates of change of the two populations at a specific point in time and space
  model <- function(t, y, parms) {
    # t : time
    # y : state vector (N1(t), N2(t))
    # parms : the parameters
    with(as.list(c(y, parms)), {
      dN1 <- N1 * (r1 - alpha11 * N1 + gamma12 * N2)
      dN2 <- N2 * (r2 + gamma21 * N1 - alpha22 * N2)
      list(c(dN1, dN2))
    })
  }

  # Defining the time span of the simulation
  times <- seq(0, 120, 0.01)

  # renaming the input parameters vector to 'parms' to use it inside the model
  parms <- parameters

  # Run simulation
  # ode takes the initial_condition, the times, the model function, and the parms
  # integrates the system of ODEs across the time points and returns a matrix
  # column 1 : time
  # column 2 : values of N_1 over time
  # column 3 : values of N_2 over time
  out <- ode(initial_condition, times, model, parms, rtol = 1e-5, atol = 1e-6)

  # -------- Calculating equilibrium points ----------------------------------

  # Trivial equilibria

  # Both species are extinct
  eq1 <- c(0, 0)
  # N1 reaches its carrying capacity when alone
  # N2 is extinct
  eq2 <- c(parms["r1"]/parms["alpha11"], 0)
  # N2 reaches its carrying capacity when alone
  # N1 extinct
  eq3 <- c(0, parms["r2"]/parms["alpha22"])

  # Coexistence equilibrium
  # To calculate the non-trivial equilibrium point, we rearrange the set of coupled ODEs as a linear system Ax = b, where :
  # A : 2x2 matrix containing the relevant parameters row1 : [alpha11, gamma12], row2:[gamma21, alpha22]
  # N : vector [N1,N2]
  # b : vector [r1, r2]
  A <- matrix(c(parms["alpha11"], parms["gamma12"],
                parms["gamma21"], parms["alpha22"]),
              nrow = 2, byrow = TRUE)
  b <- c(parms["r1"], parms["r2"])
  # Solve the linear system
  eq4 <- solve(A, b)

  # --- Jacobian and stability analysis at the coexistence equilibrium ---

  # Extract equilibrium values
  N1_star <- eq4[1]
  N2_star <- eq4[2]

  # Compute Jacobian matrix at (N1*, N2*) i.e. the equilibrium point
  # This was done to eliminate doubt about the nature of the equilibrium
  J <- matrix(c(
    parms["r1"] - 2 * parms["alpha11"] * N1_star + parms["gamma12"] * N2_star,
    parms["gamma12"] * N1_star,
    parms["gamma21"] * N2_star,
    parms["r2"] + parms["gamma21"] * N1_star - 2 * parms["alpha22"] * N2_star
  ), nrow = 2, byrow = TRUE)

  # Calculating the eigenvalues
  # all real parts negative	-> Stable node
  # in every other case -> unstable
  eig_vals <- eigen(J)$values

  print(eig_vals)


  # Determine stability
  stability_label <- if (all(Re(eig_vals) < 0)) "Stable" else "Unstable"

  # ---------------- Plotting ------------------------------------------------

  # Creating thr phase plane
  # Plotting the trajectories of N1 and N2 in the phase space
  plot(out[, "N1"], out[, "N2"],
       type = "l", col = "red", lwd = 2,
       xlab = "Species 1 (N1)", ylab = "Species 2 (N2)",
       main = paste("Lotka-Volterra Mutualistic Species\nInitial: N1 =", initial_condition["N1"],
                    ", N2 =", initial_condition["N2"], "\nEquilibrium:", stability_label),
       xlim = c(0, max(out[, "N1"]) * 1.1),
       ylim = c(0, max(out[, "N2"]) * 1.1))

  # Add isoclines
  # N1-isocline (dN1/dt = 0): N1 = 0 or r1 - alpha11*N1 + gamma12*N2 = 0
  curve((parms["r1"] - parms["alpha11"] * x)/parms["gamma12"],
        from = 0, to = max(out[, "N1"]) * 1.1,
        col = "blue", lwd = 2, lty = 2, add = TRUE)
  # N2-isocline (dN2/dt = 0): N2 = 0 or r2 + gamma21*N1 - alpha22*N2 = 0
  curve((parms["r2"] - parms["gamma21"] * x)/parms["alpha22"],
        from = 0, to = max(out[, "N1"]) * 1.1,
        col = "darkgreen", lwd = 2, lty = 2, add = TRUE)

  # Adding the equilibrium points
  points(eq2[1], eq2[2], pch = 16, cex = 1.5, col = "#8B4513")
  points(eq3[1], eq3[2], pch = 16, cex = 1.5, col = "#228B22")
  points(eq4[1], eq4[2], pch = 16, cex = 1.8, col = "purple")
  # Label the coexistence equilibrium as Stable/Unstable
  text(eq4[1], eq4[2], labels = stability_label, pos = 4, col = ifelse(stability_label == "Stable", "darkgreen", "red"))

  # Adding the vector field
  N1_seq <- seq(0, max(out[, "N1"]) * 1.1, length.out = 15)
  N2_seq <- seq(0, max(out[, "N2"]) * 1.1, length.out = 15)
  grid <- expand.grid(N1 = N1_seq, N2 = N2_seq)

  for(i in 1:nrow(grid)) {
    N1 <- grid$N1[i]
    N2 <- grid$N2[i]
    derivs <- model(0, c(N1 = N1, N2 = N2), parms)[[1]]
    # Scaling arrows optimise visualisation
    arrows(N1, N2,
           N1 + derivs[1] * scale,
           N2 + derivs[2] * scale,
           length = 0.05, col = "gray")
  }

  # Adding the legend
  legend("topright",
         legend = c("Trajectory", "N1-isocline", "N2-isocline",
                    "Equilibrium (N1 only)", "Equilibrium (N2 only)", "Equilibrium (coexistence)"),
         col = c("red", "blue", "darkgreen", "#8B4513", "#228B22", "purple"),
         lty = c(1, 2, 2, NA, NA, NA),
         lwd = c(2, 2, 2, NA, NA, NA),
         pch = c(NA, NA, NA, 16, 16, 16),
         bty = "n")
}

# (1) there exixsts a positive equilibrium points which is stable
lotka_volterra_mutalistic(
  initial_condition = c(N1 = 15, N2 = 25),
  parameters = c(
    r1 = 12, r2 = 10,
    alpha11 = 1.6, alpha22 = 1.7,
    gamma12 = 0.5, gamma21 = 0.3
  )
)

# (2)  there exists a positive equilibrium points, but which is unstable
lotka_volterra_mutalistic(
  initial_condition = c(N1 = 0.5, N2 = 0.5),
  parameters = c(
    r1 = 1.2, r2 = 1.0,
    alpha11 = 1.6, alpha22 = 1.35,
    gamma12 = 0.6, gamma21 = 0.5
  )
)

# moves through the analytical equilibrium point and away from it -> unstable

# Discussion of the results

# With higher intrinsic growth rates (r1 = 12, r2 = 10) and moderate mutualistic
# coefficients (gamma12 = 0.5, gamma21 = 0.3), the system settles at a positive coexistence
# The two isoclines intersect at a positive equilibrium indicated in purple.

# The phase plane shows the isoclines intersecting in the positive quadrant, with
# the vector field and trajectory both clearly converging towards the equilibrium point.
# The stability analysis confirmed this, as the Jacobian matrix had eigenvalues
# with negative real parts (lambda_1 = -6.768772 < 0, lambda_2 = -2.600878),
# indicating a locally stable node.

# For the totality of the system this means that :
# Both species benefit from each other but not excessively. Density-dependent
# self-limitation stabilises the system.



# In the second example (r1 = 1.2, r2 = 1.0, alpha11 = 1.6, alpha22 = 1.35,
# gamma12 = 0.6, gamma21 = 0.5), the system also has a positive coexistence equilibrium,
# but the behaviour is very different:


# Here, the mutualistic terms (gamma12 = 0.6, gamma21 = 0.5) are stronger relative
# to the intrinsic growth rates and self-limitation. Despite the trajectory starting
# close to the equilibrium point (purple), it quickly diverges. The vector field flows
# away from the intersection point, and the Jacobian matrix reveals one
# eigenvalue with a positive real part (lambda_1 = -0.50350466 < 0,
# lambda_2 = 0.09382725 > 0), indicating that this coexistence point
# is locally unstable.


# For the totality of the system this means that :
# Mutualistic interactions are too strong compared to self-limitation. This means
# that the mutualistic feedback loops can overpower self regulation and a small
# small can already be enough to push the system into uncontrolled growth also
# called runaway amplification (in this case) or collapse.

