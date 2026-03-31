# Load the required library
library(markovchain)

# --- Define the Transition Matrix ---
states_A1 <- as.character(1:5)
P_A1 <- matrix(c(
  1.0, 0,   0,   0,   0,
  0.5, 0,   0,   0,   0.5,
  0.2, 0.2, 0.2, 0.2, 0.2,
  0,   0,   1.0, 0,   0,
  0,   0,   0,   1.0, 0
), nrow = 5, byrow = TRUE)

mc_A1 <- new("markovchain", states = states_A1, byrow = TRUE, transitionMatrix = P_A1, name = "A1")

#(a) Plot and Classify 
plot(mc_A1, main = "Diagram of 5-State Markov Chain")

print("Transient Classes:")
transientClasses(mc_A1)

print("Recurrent Classes:")
recurrentClasses(mc_A1)

print("Absorbing States:")
absorbingStates(mc_A1)

print("Period of States:")
period(mc_A1)

 #(b) Simulate Three Trajectories
set.seed(42) # For reproducibility
for (i in 1:3) {
  start_state <- sample(states_A1, 1)
  sim <- rmarkovchain(n = 20, object = mc_A1, t0 = start_state)
  cat("Trajectory", i, "starting at", start_state, ":\n")
  print(sim)
}

#(c) Steady-State Probabilities
print("Steady-State Distribution:")
steadyStates(mc_A1)

print("Is the chain ergodic?")
is.irreducible(mc_A1) # If FALSE, it's not ergodic.

#(d) Plot Unconditional Probabilities vs Time
# Assume an equal starting probability for all 5 states
initial_dist <- c(0.2, 0.2, 0.2, 0.2, 0.2)
n_steps <- 20
prob_matrix <- matrix(NA, nrow = n_steps, ncol = 5)

for (n in 1:n_steps) {
  prob_matrix[n, ] <- initial_dist * (mc_A1 ^ n)
}

matplot(1:n_steps, prob_matrix, type = "l", lty = 1, col = 1:5,
        xlab = "Time (n)", ylab = "Probability", main = "Convergence to Steady State")
legend("right", legend = states_A1, col = 1:5, lty = 1)

###A2#####
states_A2 <- as.character(1:7)
P_A2 <- matrix(c(
  0,   1,   0,   0,   0,   0,   0,
  1,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0.4, 0.2, 0.2, 0.2,
  0,   0,   0,   0,   0.2, 0.4, 0.4,
  0.3, 0,   0,   0.1, 0.3, 0.1, 0.2,
  0,   0,   0,   0.2, 0.2, 0.3, 0.3,
  0,   0,   0,   0.5, 0.2, 0.2, 0.1
), nrow = 7, byrow = TRUE)

mc_A2 <- new("markovchain", states = states_A2, byrow = TRUE, transitionMatrix = P_A2, name = "A2")

# (a) Plot
plot(mc_A2, main = "Diagram of 7-State Markov Chain")

#(b) Classify States
print("Recurrent Classes:")
recurrentClasses(mc_A2)

print("Transient Classes:")
transientClasses(mc_A2)

print("Absorbing States:")
absorbingStates(mc_A2) # Empty if none exist

print("Period:")
period(mc_A2)

#(c) Simulate Two Trajectories ---
set.seed(123)
for (i in 1:2) {
  start_state <- sample(states_A2, 1)
  sim <- rmarkovchain(n = 30, object = mc_A2, t0 = start_state)
  cat("Trajectory", i, "starting at", start_state, ":\n")
  print(sim)
}

#(d) Limiting Probabilities
print("Limiting Probabilities:")
steadyStates(mc_A2)

print("Is the chain ergodic?")
is.irreducible(mc_A2)

#Define the Matrices
states_traffic <- c("light", "heavy", "jammed")

# Using your corrected value of 0.4 for row 2, column 2 
P_1to4 <- matrix(c(
  0.4, 0.4, 0.2,
  0.3, 0.4, 0.3,
  0,   0.1, 0.9
), nrow = 3, byrow = TRUE)

P_4to6 <- matrix(c(
  0.1, 0.5, 0.4,
  0.1, 0.3, 0.6,
  0,   0.1, 0.9
), nrow = 3, byrow = TRUE) 

mc_1to4 <- new("markovchain", states = states_traffic, byrow = TRUE, transitionMatrix = P_1to4)
mc_4to6 <- new("markovchain", states = states_traffic, byrow = TRUE, transitionMatrix = P_4to6)

#(a) Distribution at 6 PM
# 1 PM to 4 PM = 3 hours = 180 mins. Transitions = 180 / 20 = 9 steps.
# 4 PM to 6 PM = 2 hours = 120 mins. Transitions = 120 / 20 = 6 steps.

initial_state <- c(1, 0, 0) # Starts at "light"

# Matrix multiplication for n-steps
dist_4PM <- initial_state * (mc_1to4 ^ 9)
dist_6PM <- dist_4PM * (mc_4to6 ^ 6)

print("Analytical Distribution at 6 PM:")
print(dist_6PM)

# --- (b) Simulate 10,000 Trajectories (Corrected) ---
n_sims <- 10000
final_states <- character(n_sims)

set.seed(99)
for (i in 1:n_sims) {
  # Phase 1: 1 PM to 4 PM (9 steps)
  current_state <- "light"
  for (step in 1:9) {
    # Directly extract the row of probabilities for the current state
    probs <- mc_1to4@transitionMatrix[current_state, ]
    current_state <- sample(states_traffic, 1, prob = probs)
  }
  
  # Phase 2: 4 PM to 6 PM (6 steps)
  for (step in 1:6) {
    # Directly extract the row of probabilities from the second matrix
    probs <- mc_4to6@transitionMatrix[current_state, ]
    current_state <- sample(states_traffic, 1, prob = probs)
  }
  
  final_states[i] <- current_state
}

# Calculate proportions from simulation
simulated_dist <- prop.table(table(final_states))

print("Simulated Distribution at 6 PM:")
print(simulated_dist)

