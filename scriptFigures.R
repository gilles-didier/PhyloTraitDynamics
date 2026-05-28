library(PhyloTraitDynamics)
library(ape)
library(TreeSim)

#Simulation Brownian realisation on fixed tree Figure 1
tree <- read.tree(text = "(((A:0.2,B:0.2):0.3,C:0.4):0.5,(D:0.8, E:0.4):0.2);")
df <- fixed_tree_simulate_brownian_realization(tree, sigma2 = 100, time_step=0.001, seed=2222)
pdf(file = "simulationExt_R.pdf", width = 10, height = 7.5)
fixed_tree_plot_brownian_realization(df, cex=1.3)
dev.off()

#Theoretical distribution on fixed tree Figure 2
set.seed(1)
tr <- TreeSim::sim.bd.age(
  age = 3, lambda = 1.5, mu = 1.,
  numbsim = 1, mrca = TRUE, complete = TRUE
)[[1]]
df <- fixed_tree_compute_theoretical_summary(
  tree   = tr,
  sigma2 = 1,
  time_end = NULL,
  time_step = 0.001)
pdf(file = "SimuTreeDist_R.pdf", width = 10, height = 7.5)
fixed_tree_plot_theoretical_summary(df, cex = 1.3)
dev.off()

#MRCA age distribution Figure 3
step <- 0.005
res2_1 <- mrca_age_compute_dynamics(birth = 2,  death = 1, time_end = 5, time_step=step)
res2_2 <- mrca_age_compute_dynamics(birth = 2,  death = 2, time_end = 5, time_step=step)
res2_3 <- mrca_age_compute_dynamics(birth = 2,  death = 3, time_end = 5, time_step=step)
pdf(file = "MRCA_age_R.pdf", width = 10, height = 15)
op <- par(mfrow = c(3, 1))
par(cex = 1.3, mar = c(2, 4, 2, 1), oma = c(4, 3, 0, 0))
mrca_age_plot_dynamics(res2_1, main=expression(paste("Supercritical, ", lambda, " = 2, ", mu, " = 1")), xlab = "", ylab = "")
mrca_age_plot_dynamics(res2_2, main=expression(paste("Critical, ", lambda, " = 2, ", mu, " = 2")), xlab = "", ylab = "", add_legend = FALSE)
mrca_age_plot_dynamics(res2_3, main=expression(paste("Subcritical, ", lambda, " = 2, ", mu, " = 3")), xlab = "", ylab = "", add_legend = FALSE)
mtext("MRCA age", side = 2, outer = TRUE, line = 0., cex = 1.3)
mtext("Time", side = 1, outer = TRUE, line = 1., cex = 1.3)
par(op)
dev.off()

# Simulations computation / It may take a while
step <- 0.05
simSur <- birth_death_brownian_simulate(birth = 2, death = 1, sigma2 = 1, time_end = 5, time_step = step, B = 50000, x0 = 0, seed = 1)
simCri <- birth_death_brownian_simulate(birth = 2, death = 2, sigma2 = 1, time_end = 5, time_step = step, B = 500000, x0 = 0, seed = 1)
simSub <- birth_death_brownian_simulate(birth = 2, death = 3, sigma2 = 1, time_end = 5, time_step = step, B = 5000000, x0 = 0, seed = 1)

meanSur <- empirical_mean_compute_variance(birth = 2,  death = 1, time_end = 5, time_step=step)
meanCri <- empirical_mean_compute_variance(birth = 2,  death = 2, time_end = 5, time_step=step)
meanSub <- empirical_mean_compute_variance(birth = 2,  death = 3, time_end = 5, time_step=step)

varUnSur <- empirical_variance_compute_expectation(birth = 2,  death = 1, time_end = 5, time_step=step, method ="numeric", conditioning="none")
varUnCri <- empirical_variance_compute_expectation(birth = 2,  death = 2, time_end = 5, time_step=step, method ="numeric", conditioning="none")
varUnSub <- empirical_variance_compute_expectation(birth = 2,  death = 3, time_end = 5, time_step=step, method ="numeric", conditioning="none")

varCoSur <- empirical_variance_compute_expectation(birth = 2,  death = 1, time_end = 5, time_step=step, method ="numeric", conditioning="survival")
varCoCri <- empirical_variance_compute_expectation(birth = 2,  death = 2, time_end = 5, time_step=step, method ="numeric", conditioning="survival")
varCoSub <- empirical_variance_compute_expectation(birth = 2,  death = 3, time_end = 5, time_step=step, method ="numeric", conditioning="survival")

#Empirical mean variance, conditioned on survival, Figure 4
pdf(file = "Mean_Variance_Empirical_Mean_R.pdf", width = 10, height = 15)
op <- par(mfrow = c(3, 1))
	par(cex = 1.3, mar = c(2, 4, 2, 1), oma = c(4, 3, 0, 0))
	ylim <- range(simSur$summary$empirical_mean_variance_empirical_cond_survival, meanSur$empirical_mean_variance, na.rm = TRUE)
	plot(simSur$time, simSur$summary$empirical_mean_variance_empirical_cond_survival, col="gray75", lwd = 7, type="l", ylim = ylim, ylab = "", xlab = "", main=expression(paste("Supercritical, ", lambda, " = 2, ", mu, " = 1")))
	lines(meanSur$time, meanSur$empirical_mean_variance, col = "purple4", lty = 2, lwd = 5)
	legend("topleft", legend = c("Simulations", "Theory"), col = c("gray75","purple4"), lty = c(1,5), lwd = c(7, 5), bty = "n")
	ylim <- range(simCri$summary$empirical_mean_variance_empirical_cond_survival, meanCri$empirical_mean_variance, na.rm = TRUE)
	plot(simCri$time, simCri$summary$empirical_mean_variance_empirical_cond_survival, col="gray75", lwd = 7, type="l", ylim = ylim, ylab = "", xlab = "", main=expression(paste("Critical, ", lambda, " = 2, ", mu, " = 2")))
	lines(meanCri$time, meanCri$empirical_mean_variance, col = "purple4", lty = 2, lwd = 5)
	ylim <- range(simSub$summary$empirical_mean_variance_empirical_cond_survival, meanSub$empirical_mean_variance, na.rm = TRUE)
	plot(simSub$time, simSub$summary$empirical_mean_variance_empirical_cond_survival, col="gray75", lwd = 7, type="l", ylim = ylim, ylab = "", xlab = "", main=expression(paste("Subcritical, ", lambda, " = 2, ", mu, " = 3")))
	lines(meanSub$time, meanSub$empirical_mean_variance, col = "purple4", lty = 2, lwd = 5)
mtext("Empirical mean\nvariance", side = 2, outer = TRUE, line = 0., cex = 1.3)
mtext("Time", side = 1, outer = TRUE, line = 1., cex = 1.3)
dev.off()


#Expected empirical variance, unconditioned, Figure 5
pdf(file = "Mean_Expected_Empirical_Variance_Un_R.pdf", width = 10, height = 15)
op <- par(mfrow = c(3, 1))
	par(cex = 1.3, mar = c(2, 4, 2, 1), oma = c(4, 3, 0, 0))
	ylim <- range(simSur$summary$empirical_variance_expectation_empirical, varUnSur$empirical_variance_expectation, na.rm = TRUE)
	plot(simSur$time, simSur$summary$empirical_variance_expectation_empirical, col="gray75", lwd = 7, type="l", ylim = ylim, ylab = "", xlab = "", main=expression(paste("Supercritical, ", lambda, " = 2, ", mu, " = 1")))
	lines(varUnSur$time, varUnSur$empirical_variance_expectation, col = "darkorange3", lty = 2, lwd = 5)
	legend("topleft", legend = c("Simulations", "Theory"), col = c("gray75","darkorange3"), lty = c(1,5), lwd = c(7, 5), bty = "n")
	ylim <- range(simCri$summary$empirical_variance_expectation_empirical, varUnCri$empirical_variance_expectation, na.rm = TRUE)
	plot(simCri$time, simCri$summary$empirical_variance_expectation_empirical, col="gray75", lwd = 7, type="l", ylim = ylim, ylab = "", xlab = "", main=expression(paste("Critical, ", lambda, " = 2, ", mu, " = 2")))
	lines(varUnCri$time, varUnCri$empirical_variance_expectation, col = "darkorange3", lty = 2, lwd = 5)
	ylim <- range(simSub$summary$empirical_variance_expectation_empirical, varUnSub$empirical_variance_expectation, na.rm = TRUE)
	plot(simSub$time, simSub$summary$empirical_variance_expectation_empirical, col="gray75", lwd = 7, type="l", ylim = ylim, ylab = "", xlab = "", main=expression(paste("Subcritical, ", lambda, " = 2, ", mu, " = 3")))
	lines(varUnSub$time, varUnSub$empirical_variance_expectation, col = "darkorange3", lty = 2, lwd = 5)
mtext("Expected\nempirical variance", side = 2, outer = TRUE, line = 0., cex = 1.3)
mtext("Time", side = 1, outer = TRUE, line = 1., cex = 1.3)
dev.off()

#Expected empirical variance, conditioned on survival, Figure 6
pdf(file = "Mean_Expected_Empirical_Variance_Co_R.pdf", width = 10, height = 15)
op <- par(mfrow = c(3, 1))
	par(cex = 1.3, mar = c(2, 4, 2, 1), oma = c(4, 3, 0, 0))
	ylim <- range(simSur$summary$empirical_variance_expectation_empirical_cond_survival, varCoSur$empirical_variance_expectation, na.rm = TRUE)
	plot(simSur$time, simSur$summary$empirical_variance_expectation_empirical_cond_survival, col="gray75", lwd = 7, type="l", ylim = ylim, ylab = "", xlab = "", main=expression(paste("Supercritical, ", lambda, " = 2, ", mu, " = 1")))
	lines(varCoSur$time, varCoSur$empirical_variance_expectation, col = "darkorange3", lty = 2, lwd = 5)
	legend("topleft", legend = c("Simulations", "Theory"), col = c("gray75","darkorange3"), lty = c(1,5), lwd = c(7, 5), bty = "n")
	ylim <- range(simCri$summary$empirical_variance_expectation_empirical_cond_survival, varCoCri$empirical_variance_expectation, na.rm = TRUE)
	plot(simCri$time, simCri$summary$empirical_variance_expectation_empirical_cond_survival, col="gray75", lwd = 7, type="l", ylim = ylim, ylab = "", xlab = "", main=expression(paste("Critical, ", lambda, " = 2, ", mu, " = 2")))
	lines(varCoCri$time, varCoCri$empirical_variance_expectation, col = "darkorange3", lty = 2, lwd = 5)
	ylim <- range(simSub$summary$empirical_variance_expectation_empirical_cond_survival, varCoSub$empirical_variance_expectation, na.rm = TRUE)
	plot(simSub$time, simSub$summary$empirical_variance_expectation_empirical_cond_survival, col="gray75", lwd = 7, type="l", ylim = ylim, ylab = "", xlab = "", main=expression(paste("Subcritical, ", lambda, " = 2, ", mu, " = 3")))
	lines(varCoSub$time, varCoSub$empirical_variance_expectation, col = "darkorange3", lty = 2, lwd = 5)
mtext("Expected\nempirical variance", side = 2, outer = TRUE, line = 0., cex = 1.3)
mtext("Time", side = 1, outer = TRUE, line = 1., cex = 1.3)
dev.off()


#Variable rates Figure 7
compute_band_limits <- function(tmin, tmax, n_bands = 10) {
  width_raw <- (tmax - tmin) / n_bands
  step <- 10 ^ floor(log10(width_raw))
  if(step == width_raw) {
	  width = width_raw
  } else {
	  width = step*10
	  min = abs(width_raw-width)
	  if(abs(width_raw-width*0.5)<min) {
		  width = width*0.5
		  min = abs(width_raw-width)
	 }
	  if(abs(width_raw-width*0.5)<min) {
		  width = width*0.5
		  min = abs(width_raw-width)
	 }
	  if(abs(width_raw-step)<min) {
		  width = step
		  min = abs(width_raw-width)
	 }
  }
  first_start <- floor(tmin / width) * width
  starts <- seq(first_start, tmax, by = width)
  ends <- starts + width / 2
  ends[length(ends)] <- min(ends[length(ends)], tmax)
  data.frame(start = starts, end = ends)
}

add_time_bands <- function(bands, col = grDevices::adjustcolor("grey", alpha.f = 0.2)) {
  usr <- graphics::par("usr")
  for (i in seq_len(nrow(bands))) {
    graphics::rect(bands$start[i], usr[3], bands$end[i], usr[4], col = col, border = NA)
  }
}
fun_bump <- function(a, b, center = 2.5, width = 0.6) {
  stopifnot(a >= 0, b >= a, width > 0)
  
  function(t) {
    a + (b - a) * exp(-((t - center) / width)^2)
  }
}
	birth <- fun_bump(0.2, 10, width = 0.1)
	death <- function(t) { 
		rep(0., length(t))}
	time_end = 5
	step = 0.05
	time_band_col = grDevices::adjustcolor("blue", alpha.f = 0.05)
	mrca <- mrca_age_compute_dynamics(birth,  death, time_end = 5, time_step=step)
	mean <- empirical_mean_compute_variance(birth,  death, time_end = 5, time_step=step)
	var <- empirical_variance_compute_expectation(birth,  death, time_end = 5, time_step=step, method ="numeric", conditioning="survival")
	nb_simul = 1000
	simVar <- birth_death_brownian_simulate(birth = birth, death = death, sigma2 = 1, time_end = time_end, time_step = step, B = nb_simul, x0 = 0, seed = 1)
	bands <- compute_band_limits(0, time_end)
	
	pdf(file = "VariableRates_Paths_R.pdf", width = 10, height = 10)
	layout(matrix(1:4, ncol = 1), heights = c(0.8, 1, 1, 1.5))
	par(cex = 1.3)
	par(mar = c(0, 6, 2, 1))
	plot(birth, 0, 5, col="forestgreen", type="l", axes = FALSE, ylim = c(-1, 11), ylab = "Diversification\nrates", n = 1000)
	plot(death, 0, 5, col= "brown", type="l", add =TRUE)
	legend("topleft", legend = c("Birth", "Death"), col = c("forestgreen","brown"), lty = c(1,1), lwd = c(1, 1), bty = "n")
	Axis(side=2)
	add_time_bands(bands, time_band_col)
	par(mar = c(0, 6, 0, 1))
	mrca_age_plot_dynamics(mrca, axes = FALSE, ylab = "MRCA\nage")
	Axis(side=2)
	add_time_bands(bands, time_band_col)
	plot(simVar$time, simVar$summary$empirical_mean_variance_empirical_cond_survival, col = "gray75", type='l', axes = FALSE,  lwd = 7, ylab = "Empirical mean\nvariance", xlab = "")
	lines(mean$time, mean$empirical_mean_variance, type='l', col = "purple4", lwd = 3, lty = 2)
	legend("topleft", legend = c("Simulations", "Theory"), col = c("gray75", "purple4"), lty = c(1,5), lwd = c(7, 5), bty = "n")
	ylim3 <- range(c(mean$empirical_mean_variance,  simVar$summary$empirical_mean_variance_empirical_cond_survival), na.rm = TRUE)
	yticks <- pretty(ylim3)
	yticks <- yticks[yticks >= ylim3[1] & yticks <= ylim3[2]]
	Axis(side=2,     at = yticks,
    labels = c(yticks[-length(yticks)], ""), las = 1)
	add_time_bands(bands, time_band_col)
	par(mar = c(6, 6, 0, 1))
	plot(simVar$time, simVar$summary$empirical_variance_expectation_empirical_cond_survival, col = "gray75", type='l', axes = FALSE,  lwd = 7, ylab = "Expected\nempirical variance", xlab = "Time")
	lines(var$time, var$empirical_variance_expectation, col = "darkorange3", lwd = 3, lty = 2)
	legend("topleft", c("Simulations", "Theory"), col = c("gray75", "darkorange3"), lty = c(1,5), lwd = c(7, 5), bty = "n")
	Axis(side=1)
	Axis(side=2)
	add_time_bands(bands, time_band_col)
dev.off()

