require(cowplot)
require(dplyr)
require(purrr)
require(xtable)
require(stringr)
require(truncnorm)
require(RColorBrewer)
require(reshape2)
require(gridGraphics)
require(plotly)
library(listviewer)
load("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Environments/bootmain.RData")
alp <- 0.05
ntau <- 11
vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
##########################
#Lower Bound
##########################
#For output elasticities that vary over percentiles of inputs
k3dLB <- array(0, c(ntau, ntau))
l3dLB <- k3dLB; m3dLB <- k3dLB
#For output elasticities that vary over productivity
kwq3dLB <- k3dLB; lwq3dLB <- k3dLB; mwq3dLB <- k3dLB
#For output elasticities that vary over percentiles of other inputs
klq3dLB <- k3dLB; kmq3dLB <- k3dLB
lkq3dLB <- k3dLB; lmq3dLB <- k3dLB
mkq3dLB <- k3dLB; mlq3dLB <- k3dLB
#For input demand functions varying over percentiles of productivity
iw3dLB <- k3dLB; lw3dLB <- k3dLB; mw3dLB <- k3dLB
#For input demand functions varying over percentiles of capital
iwkq3dLB <- k3dLB; lwkq3dLB <- k3dLB; mwkq3dLB <- k3dLB
#For non-Hicks neutral effects varying over percentiles of inputs
hk3dLB <- k3dLB
hl3dLB <- hk3dLB; hm3dLB <- hk3dLB
#For non-Hicks neutral effects that vary over percentiles of other inputs
hklq3dLB <- k3dLB; hkmq3dLB <- k3dLB
hlkq3dLB <- k3dLB; hlmq3dLB <- k3dLB
hmkq3dLB <- k3dLB; hmlq3dLB <- k3dLB
#For scale and ratio effects
hsLB <- k3dLB; hklLB <- k3dLB; hkmLB <- k3dLB; hlmLB <- k3dLB
#Persistence
omgLB <- k3dLB
#Dispersion, Skewness, and Kurtosis
dispLB <- array(0, c(ntau)); skewLB <- dispLB; kurtLB <- dispLB
##########################
#Upper Bound
##########################
#For output elasticities that vary over percentiles of inputs
k3dUB <- array(0, c(ntau, ntau))
l3dUB <- k3dUB; m3dUB <- k3dUB
#For output elasticities that vary over productivity
kwq3dUB <- k3dUB; lwq3dUB <- k3dUB; mwq3dUB <- k3dUB
#For output elasticities that vary over percentiles of other inputs
klq3dUB <- k3dUB; kmq3dUB <- k3dUB
lkq3dUB <- k3dUB; lmq3dUB <- k3dUB
mkq3dUB <- k3dUB; mlq3dUB <- k3dUB
#For input demand functions varying over percentiles of productivity
iw3dUB <- k3dUB; lw3dUB <- k3dUB; mw3dUB <- k3dUB
#For input demand functions varying over percentiles of capital
iwkq3dUB <- k3dUB; lwkq3dUB <- k3dUB; mwkq3dUB <- k3dUB
#For non-Hicks neutral effects varying over percentiles of inputs
hk3dUB <- k3dUB
hl3dUB <- hk3dUB; hm3dUB <- hk3dUB
#For non-Hicks neutral effects that vary over percentiles of other inputs
hklq3dUB <- k3dUB; hkmq3dUB <- k3dUB
hlkq3dUB <- k3dUB; hlmq3dUB <- k3dUB
hmkq3dUB <- k3dUB; hmlq3dUB <- k3dUB
#For scale and ratio effects
hsUB <- k3dUB; hklUB <- k3dUB; hkmUB <- k3dUB; hlmUB <- k3dUB
#Persistence
omgUB <- k3dUB
#Dispersion, Skewness, and Kurtosis
dispUB <- array(0, c(ntau)); skewUB <- dispUB; kurtUB <- dispUB
######################
#Confidence Bands
#####################
for (q1 in 1:ntau){
	#Dispersion
	dispLB[q1] <- as.numeric(quantile(bootresults$dispboot[q1,], probs=alp))
	dispUB[q1] <- as.numeric(quantile(bootresults$dispboot[q1,], probs=1-alp))
	#Skewness
	skewLB[q1] <- as.numeric(quantile(bootresults$skewboot[q1,], probs=alp))
	skewUB[q1] <- as.numeric(quantile(bootresults$skewboot[q1,], probs=1-alp))
	#Kurtosis
	kurtLB[q1] <- as.numeric(quantile(bootresults$kurtboot[q1,], probs=alp))
	kurtUB[q1] <- as.numeric(quantile(bootresults$kurtboot[q1,], probs=1-alp))
	for (q2 in 1:ntau){
		##################################
		#Elasticities
		##################################
		#Capital Elasticity (over capital)
		k3dLB[q1,q2] <- as.numeric(quantile(bootresults$k3dboot[q1,,][q2,], probs=alp))
		k3dUB[q1,q2] <- as.numeric(quantile(bootresults$k3dboot[q1,,][q2,], probs=1-alp))
		#Capital Elasticity (over productivity)
		kwq3dLB[q1,q2] <- as.numeric(quantile(bootresults$kwq3dboot[q1,,][q2,], probs=alp))
		kwq3dUB[q1,q2] <- as.numeric(quantile(bootresults$kwq3dboot[q1,,][q2,], probs=1-alp))
		#Capital Elasticity (over labor)
		klq3dLB[q1,q2] <- as.numeric(quantile(bootresults$klq3dboot[q1,,][q2,], probs=alp))
		klq3dUB[q1,q2] <- as.numeric(quantile(bootresults$klq3dboot[q1,,][q2,], probs=1-alp))
		#Capital Elasticity (over materials)
		kmq3dLB[q1,q2] <- as.numeric(quantile(bootresults$kmq3dboot[q1,,][q2,], probs=alp))
		kmq3dUB[q1,q2] <- as.numeric(quantile(bootresults$kmq3dboot[q1,,][q2,], probs=1-alp))
		#Labor Elasticity (over labor)
		l3dLB[q1,q2] <- as.numeric(quantile(bootresults$l3dboot[q1,,][q2,], probs=alp))
		l3dUB[q1,q2] <- as.numeric(quantile(bootresults$l3dboot[q1,,][q2,], probs=1-alp))
		#Labor Elasticity (over productivity)
		lwq3dLB[q1,q2] <- as.numeric(quantile(bootresults$lwq3dboot[q1,,][q2,], probs=alp))
		lwq3dUB[q1,q2] <- as.numeric(quantile(bootresults$lwq3dboot[q1,,][q2,], probs=1-alp))
		#Labor Elasticity (over capital)
		lkq3dLB[q1,q2] <- as.numeric(quantile(bootresults$lkq3dboot[q1,,][q2,], probs=alp))
		lkq3dUB[q1,q2] <- as.numeric(quantile(bootresults$lkq3dboot[q1,,][q2,], probs=1-alp))
		#Labor Elasticity (over materials)
		lmq3dLB[q1,q2] <- as.numeric(quantile(bootresults$lmq3dboot[q1,,][q2,], probs=alp))
		lmq3dUB[q1,q2] <- as.numeric(quantile(bootresults$lmq3dboot[q1,,][q2,], probs=1-alp))
		#Materials Elasticity (over materials)
		m3dLB[q1,q2] <- as.numeric(quantile(bootresults$m3dboot[q1,,][q2,], probs=alp))
		m3dUB[q1,q2] <- as.numeric(quantile(bootresults$m3dboot[q1,,][q2,], probs=1-alp))
		#Labor Elasticity (over productivity)
		mwq3dLB[q1,q2] <- as.numeric(quantile(bootresults$mwq3dboot[q1,,][q2,], probs=alp))
		mwq3dUB[q1,q2] <- as.numeric(quantile(bootresults$mwq3dboot[q1,,][q2,], probs=1-alp))
		#Materials Elasticity (over capital)
		mkq3dLB[q1,q2] <- as.numeric(quantile(bootresults$mkq3dboot[q1,,][q2,], probs=alp))
		mkq3dUB[q1,q2] <- as.numeric(quantile(bootresults$mkq3dboot[q1,,][q2,], probs=1-alp))
		#Materials Elasticity (over labor)
		mlq3dLB[q1,q2] <- as.numeric(quantile(bootresults$mlq3dboot[q1,,][q2,], probs=alp))
		mlq3dUB[q1,q2] <- as.numeric(quantile(bootresults$mlq3dboot[q1,,][q2,], probs=1-alp))
		##################################
		#Non-Hicks Neutral Effects
		##################################
		#Capital (over capital)
		hk3dLB[q1,q2] <- as.numeric(quantile(bootresults$hk3dboot[q1,,][q2,], probs=alp))
		hk3dUB[q1,q2] <- as.numeric(quantile(bootresults$hk3dboot[q1,,][q2,], probs=1-alp))
		#Capital (over labor)
		hklq3dLB[q1,q2] <- as.numeric(quantile(bootresults$hklq3dboot[q1,,][q2,], probs=alp))
		hklq3dUB[q1,q2] <- as.numeric(quantile(bootresults$hklq3dboot[q1,,][q2,], probs=1-alp))
		#Capital Elasticity (over materials)
		hkmq3dLB[q1,q2] <- as.numeric(quantile(bootresults$hkmq3dboot[q1,,][q2,], probs=alp))
		hkmq3dUB[q1,q2] <- as.numeric(quantile(bootresults$hkmq3dboot[q1,,][q2,], probs=1-alp))
		#Labor (over labor)
		hl3dLB[q1,q2] <- as.numeric(quantile(bootresults$hl3dboot[q1,,][q2,], probs=alp))
		hl3dUB[q1,q2] <- as.numeric(quantile(bootresults$hl3dboot[q1,,][q2,], probs=1-alp))
		#Labor (over capital)
		hlkq3dLB[q1,q2] <- as.numeric(quantile(bootresults$hlkq3dboot[q1,,][q2,], probs=alp))
		hlkq3dUB[q1,q2] <- as.numeric(quantile(bootresults$hlkq3dboot[q1,,][q2,], probs=1-alp))
		#Labor (over materials)
		hlmq3dLB[q1,q2] <- as.numeric(quantile(bootresults$hlmq3dboot[q1,,][q2,], probs=alp))
		hlmq3dUB[q1,q2] <- as.numeric(quantile(bootresults$hlmq3dboot[q1,,][q2,], probs=1-alp))
		#Materials (over materials)
		hm3dLB[q1,q2] <- as.numeric(quantile(bootresults$hm3dboot[q1,,][q2,], probs=alp))
		hm3dUB[q1,q2] <- as.numeric(quantile(bootresults$hm3dboot[q1,,][q2,], probs=1-alp))
		#Materials (over capital)
		hmkq3dLB[q1,q2] <- as.numeric(quantile(bootresults$hmkq3dboot[q1,,][q2,], probs=alp))
		hmkq3dUB[q1,q2] <- as.numeric(quantile(bootresults$hmkq3dboot[q1,,][q2,], probs=1-alp))
		#Materials (over labor)
		hmlq3dLB[q1,q2] <- as.numeric(quantile(bootresults$hmlq3dboot[q1,,][q2,], probs=alp))
		hmlq3dUB[q1,q2] <- as.numeric(quantile(bootresults$hmlq3dboot[q1,,][q2,], probs=1-alp))
		#Scale Effect
		hsLB[q1,q2] <- as.numeric(quantile(bootresults$hsboot[q1,,][q2,], probs=alp))
		hsUB[q1,q2] <- as.numeric(quantile(bootresults$hsboot[q1,,][q2,], probs=1-alp))
		#Capital/Labor Effect
		hklLB[q1,q2] <- as.numeric(quantile(bootresults$hklboot[q1,,][q2,], probs=alp))
		hklUB[q1,q2] <- as.numeric(quantile(bootresults$hklboot[q1,,][q2,], probs=1-alp))
		#Capital/Materials Effect
		hkmLB[q1,q2] <- as.numeric(quantile(bootresults$hkmboot[q1,,][q2,], probs=alp))
		hkmUB[q1,q2] <- as.numeric(quantile(bootresults$hkmboot[q1,,][q2,], probs=1-alp))
		#Labor/Materials Effect
		hlmLB[q1,q2] <- as.numeric(quantile(bootresults$hlmboot[q1,,][q2,], probs=alp))
		hlmUB[q1,q2] <- as.numeric(quantile(bootresults$hlmboot[q1,,][q2,], probs=1-alp))
		##################################
		#Input Reponse to Productivity
		##################################
		#Investment (over productivity)
		iw3dLB[q1,q2] <- as.numeric(quantile(bootresults$iw3dboot[q1,,][q2,], probs=alp))
		iw3dUB[q1,q2] <- as.numeric(quantile(bootresults$iw3dboot[q1,,][q2,], probs=1-alp))
		#Investment (over capital)
		iwkq3dLB[q1,q2] <- as.numeric(quantile(bootresults$iwkq3dboot[q1,,][q2,], probs=alp))
		iwkq3dUB[q1,q2] <- as.numeric(quantile(bootresults$iwkq3dboot[q1,,][q2,], probs=1-alp))
		#Labor (over productivity)
		lw3dLB[q1,q2] <- as.numeric(quantile(bootresults$lw3dboot[q1,,][q2,], probs=alp))
		lw3dUB[q1,q2] <- as.numeric(quantile(bootresults$lw3dboot[q1,,][q2,], probs=1-alp))
		#Labor (over capital)
		lwkq3dLB[q1,q2] <- as.numeric(quantile(bootresults$lwkq3dboot[q1,,][q2,], probs=alp))
		lwkq3dUB[q1,q2] <- as.numeric(quantile(bootresults$lwkq3dboot[q1,,][q2,], probs=1-alp))
		#Materials (over productivity)
		mw3dLB[q1,q2] <- as.numeric(quantile(bootresults$mw3dboot[q1,,][q2,], probs=alp))
		mw3dUB[q1,q2] <- as.numeric(quantile(bootresults$mw3dboot[q1,,][q2,], probs=1-alp))
		#Materials (over capital)
		mwkq3dLB[q1,q2] <- as.numeric(quantile(bootresults$mwkq3dboot[q1,,][q2,], probs=alp))
		mwkq3dUB[q1,q2] <- as.numeric(quantile(bootresults$mwkq3dboot[q1,,][q2,], probs=1-alp))
		#Persistence
		omgLB[q1,q2] <- as.numeric(quantile(bootresults$omgboot[q1,,][q2,], probs=alp))
		omgUB[q1,q2] <- as.numeric(quantile(bootresults$omgboot[q1,,][q2,], probs=1-alp))

	}
}
#############################
#Capital Elasticities
#############################
#Capital Elasticity (over capital)
kboot <- plot_ly(x=vectau, y=vectau, z=k3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=k3dUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Capital Elasticity")))
# kboot
#Capital Elasticity (over productivity)
kwqboot <- plot_ly(x=vectau, y=vectau, z=kwq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=kwq3dUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Capital Elasticity")))
# kwqboot
#Capital Elasticity (over labor)
klqboot <- plot_ly(x=vectau, y=vectau, z=klq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=klq3dUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Capital Elasticity")))
# klqboot
#Capital Elasticity (over materials)
kmqboot <- plot_ly(x=vectau, y=vectau, z=kmq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=kmq3dUB, colorscale="YlOrRd") %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Capital Elasticity")))
# kmqboot
#Capital Non-Hicks Neutral Effect (over capital)
hkboot <- plot_ly(x=vectau, y=vectau, z=hk3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=hk3dUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Capital Efficiency")))
# hkboot
#Capital Non-Hicks Neutral Effect (over labor)
hklqboot <- plot_ly(x=vectau, y=vectau, z=hklq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=hklq3dUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Capital Efficiency")))
# hklqboot
#Capital Non-Hicks Neutral Effect (over materials)
hkmqboot <- plot_ly(x=vectau, y=vectau, z=hkmq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=hkmq3dUB, colorscale="YlOrRd") %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Capital Efficiency")))
# kmqboot
#############################
#Labor Elasticities
#############################
#Labor Elasticity (over labor)
lboot <- plot_ly(x=vectau, y=vectau, z=l3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=l3dUB, colorscale="YlOrRd") %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Labor Elasticity")))
# lboot
#Labor Elasticity (over productivity)
lwqboot <- plot_ly(x=vectau, y=vectau, z=lwq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=lwq3dUB, colorscale="YlOrRd") %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Labor Elasticity")))
# lwqboot
#Labor Elasticity (over capital)
lkqboot <- plot_ly(x=vectau, y=vectau, z=lkq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=lkq3dUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Labor Elasticity")))
# lkqboot
#Labor Elasticity (over materials)
lmqboot <- plot_ly(x=vectau, y=vectau, z=lmq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=lmq3dUB, colorscale="YlOrRd") %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Labor Elasticity")))
# lmqboot
#Labor Non-Hicks Neutral Effect (over labor)
hlboot <- plot_ly(x=vectau, y=vectau, z=hl3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=hl3dUB, colorscale="YlOrRd") %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Labor Efficiency")))
# hlboot
#Labor Non-Hicks Neutral Effect (over capital)
hlkqboot <- plot_ly(x=vectau, y=vectau, z=hlkq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=hlkq3dUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Labor Efficiency")))
# hlkqboot
#Labor Non-Hicks Neutral Effect (over materials)
hlmqboot <- plot_ly(x=vectau, y=vectau, z=hlmq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=hlmq3dUB, colorscale="YlOrRd") %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Labor Efficiency")))
# hlmqboot
#############################
#Materials Elasticities
#############################
#Materials Elasticity (over materials)
mboot <- plot_ly(x=vectau, y=vectau, z=m3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=m3dUB, colorscale="YlOrRd") %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Materials Elasticity")))
# mboot
#Materials Elasticity (over productivity)
mwqboot <- plot_ly(x=vectau, y=vectau, z=mwq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=mwq3dUB, colorscale="YlOrRd") %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Materials Elasticity")))
# mwqboot
#Materials Elasticity (over capital)
mkqboot <- plot_ly(x=vectau, y=vectau, z=mkq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=mkq3dUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Materials Elasticity")))
# mkqboot
#Materials Elasticity (over labor)
mlqboot <- plot_ly(x=vectau, y=vectau, z=mlq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=mlq3dUB, colorscale="YlOrRd") %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Materials Elasticity")))
# mlqboot
#Materials Non-Hicks Neutral Effect (over materials)
hmboot <- plot_ly(x=vectau, y=vectau, z=hm3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=hm3dUB, colorscale="YlOrRd") %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Materials Efficiency")))
# hmboot
#Materials Non-Hicks Neutral Effect (over capital)
hmkqboot <- plot_ly(x=vectau, y=vectau, z=hmkq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=hmkq3dUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Materials Efficiency")))
# hmkqboot
#Materials Non-Hicks Neutral Effect (over labor)
hmlqboot <- plot_ly(x=vectau, y=vectau, z=hmlq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=hmlq3dUB, colorscale="YlOrRd") %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Materials Efficiency")))
# hmlqboot
#############################
#Scale and Ratio Effects
#############################
#Scale Effect
hsboot <- plot_ly(x=vectau, y=vectau, z=hsLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-scale: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=hsUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output", titlefont=list(size=30), tickfont=list(size=16)), yaxis=list(title="𝛕-scale", titlefont=list(size=30), tickfont=list(size=16)), zaxis=list(title="Scale Effect", titlefont=list(size=30), tickfont=list(size=16)))) 
# hsboot
#Capital/Labor Effect
hklboot <- plot_ly(x=vectau, y=vectau, z=hklLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-(k-l): %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=hklUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-(k-l)"), zaxis=list(title="Capital/Labor Effect"))) 
# hklboot
#Capital/Materials Effect
hkmboot <- plot_ly(x=vectau, y=vectau, z=hkmLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-(k-m): %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=hkmUB, colorscale="YlOrRd") %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-(k-m)"), zaxis=list(title="Capital/Materials Effect"))) 
# hkmboot
#Labor/Materials Effect
hlmboot <- plot_ly(x=vectau, y=vectau, z=hlmLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-(l-m): %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=hlmUB, colorscale="YlOrRd") %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-(l-m)"), zaxis=list(title="Labor/Materials Effect"))) 
# hlmboot
#######################################
#Input Demand Response to Productivity
#######################################
#Investment (over productivity)
iwboot <- plot_ly(x=vectau, y=vectau, z=iw3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-investment<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=iw3dUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-investment"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Investment Productivity"))) 
# iwboot
#Investment (over capital)
iwkqboot <- plot_ly(x=vectau, y=vectau, z=iwkq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-investment<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=iwkq3dUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-investment"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Investment Productivity"))) 
# iwkqboot
#Labor (over productivity)
lwboot <- plot_ly(x=vectau, y=vectau, z=lw3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-labor<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=lw3dUB, colorscale="YlOrRd") %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-labor"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Labor Productivity"))) 
# lwboot
#Labor (over capital)
lwkqboot <- plot_ly(x=vectau, y=vectau, z=lwkq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-labor<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=lwkq3dUB, colorscale="YlOrRd") %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-labor"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Labor Productivity"))) 
# lwkqboot
#Materials (over productivity)
mwboot <- plot_ly(x=vectau, y=vectau, z=mw3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-materials<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=mw3dUB, colorscale="YlOrRd") %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-materials"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Materials Productivity"))) 
# mwboot
#Materials (over capital)
mwkqboot <- plot_ly(x=vectau, y=vectau, z=mwkq3dLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-materials<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=mwkq3dUB, colorscale="YlOrRd") %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-materials"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Materials Productivity"))) 
# mwkqboot
#Combined
lmiwboot <- subplot(iwboot, lwboot, mwboot) %>% layout(scene=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-investment", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Investment", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
lmiwkqboot <- subplot(iwkqboot, lwkqboot, mwkqboot) %>% layout(scene=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-investment", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Investment", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
#Annotations
latexannotationsw <- list(list(x=0.08, y=0.75, text="(a) Investment Response", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Response", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.92, y=0.75, text="(c) Materials Response", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationsw <- list(list(x=0.07, y=0.75, text="(a) Investment Response", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Response", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.92, y=0.75, text="(c) Materials Response", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Add
lmiwplotboot_latex <- lmiwboot %>% layout(annotations=latexannotationsw)
lmiwkqplotboot_latex <- lmiwkqboot %>% layout(annotations=latexannotationsw)
#Plot
lmiwplotboot_latex
lmiwkqplotboot_latex
#######################################
#Productivity
#######################################
#Persistence
omgboot <- plot_ly(x=vectau, y=vectau, z=omgLB, colorscale="YlGnBu", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-innovation<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% add_surface(z=omgUB, colorscale="YlOrRd") %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-innovation", titlefont=list(size=20), tickfont=list(size=16)), yaxis=list(title="𝛕-productivity", titlefont=list(size=20), tickfont=list(size=16)), zaxis=list(title="Persistence", titlefont=list(size=20), tickfont=list(size=16)))) 
# omgboot
#Dispersion
dispboot <- plot_ly(x=vectau, y=dispLB, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black"), name=" ", hovertemplate = paste("<i>𝛕-productivity<i>: %{x:.2f}", "<br>Estimate: %{y:.3f}")) %>% add_trace(y=dispUB, line=list(color="black"), fill="tonexty", fillcolor="rgba(100,100,100,0.3)") %>% layout(xaxis=list(title="𝛕-productivity", titlefont=list(size=30), tickfont=list(size=25)), yaxis=list(title="Conditional Dispersion", titlefont=list(size=18), tickfont=list(size=25), range=list(min(dispLB)-0.1*min(dispLB), max(dispUB)+0.1*max(dispUB))))
# dispboot 
#Skewness
skewboot <- plot_ly(x=vectau, y=skewLB, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black"), name=" ", hovertemplate = paste("<i>𝛕-productivity<i>: %{x:.2f}", "<br>Estimate: %{y:.3f}")) %>% add_trace(y=skewUB, line=list(color="black"), fill="tonexty", fillcolor="rgba(100,100,100,0.3)")%>% layout(xaxis=list(title="𝛕-productivity", titlefont=list(size=30), tickfont=list(size=25)), yaxis=list(title="Conditional Skewness", titlefont=list(size=18), tickfont=list(size=25), range=list(min(skewLB)-0.5*abs(min(skewLB)), max(skewUB)+0.5*abs(max(skewUB)))))
# skewboot
#Kurtosis
kurtboot <- plot_ly(x=vectau, y=kurtLB, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black"), name=" ", hovertemplate = paste("<i>𝛕-productivity<i>: %{x:.2f}", "<br>Estimate: %{y:.3f}")) %>% add_trace(y=kurtUB, line=list(color="black"), fill="tonexty", fillcolor="rgba(100,100,100,0.3)")%>% layout(xaxis=list(title="𝛕-productivity", titlefont=list(size=30), tickfont=list(size=25)), yaxis=list(title="Conditional Kurtosis", titlefont=list(size=18), tickfont=list(size=25), range=list(min(kurtLB)-0.1*min(kurtLB), max(kurtUB)+0.1*max(kurtUB))))
# kurtboot
omgbootannotate <- list(list(x=0.03, y=0.95, text="(a) Conditional Dispersion", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.95, text="(b) Conditional Skewness", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.95, y=0.95, text="(c) Conditional Kurtosis", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
omgdistboot <- subplot(dispboot, skewboot, kurtboot, shareX=TRUE) %>% layout(annotations=omgbootannotate)
omgdistboot


###################################
#Combined Plot.ly Elasticities
###################################
#Elasticities (over percentiles of inputs)
latexannotationsklm <- list(list(x=0.09, y=0.75, text="(a) Capital Elasticity", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Elasticity", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.95, y=0.75, text="(c) Materials Elasticity", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationsklm <- list(list(x=0.1, y=0.75, text="(a) Capital Elasticity", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Elasticity", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.9, y=0.75, text="(c) Materials Elasticity", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))	
klmboot <- subplot(kboot, lboot, mboot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
klmboot_latex <- klmboot %>% layout(annotations=latexannotationsklm)
klmboot_latex 
# Elasticities (over percentiles of productivity)
klmwboot <- subplot(kwqboot, lwqboot, mwqboot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
klmwboot_latex <- klmwboot %>% layout(annotations=latexannotationsklm)
klmwboot_latex
#Elasticities (over percentiles of other inputs)
latexannotationsab <- list(list(x=0.25, y=0.8, text="<b>(a)<b>", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.75, y=0.8, text="<b>(b)<b>", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationsab <- list(list(x=0.25, y=0.8, text="(a)", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.75, y=0.8, text="(b)", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Capital
kinpboot <- subplot(klqboot, kmqboot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))))
kinpboot_latex <- kinpboot %>% layout(annotations=latexannotationsab)
kinpboot_latex
#Labor
linpboot <- subplot(lkqboot, lmqboot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))))
linpboot_latex <- linpboot %>% layout(annotations=latexannotationsab)
linpboot_latex
#Materials
minpboot <- subplot(mkqboot, mlqboot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
minpboot_latex <- minpboot %>% layout(annotations=latexannotationsab)
minpboot_latex
# Non-Hicks Neutral Effects
#Annotations
latexannotationshklm <- list(list(x=0.09, y=0.75, text="(a) Capital Efficiency", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Efficiency", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.95, y=0.75, text="(c) Materials Efficiency", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationshklm <- list(list(x=0.1, y=0.75, text="(a) Capital Efficiency", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Efficiency", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.9, y=0.75, text="(c) Materials Efficiency", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Over Inputs
hklmboot <- subplot(hkboot, hlboot, hmboot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
hklmboot_latex <- hklmboot %>% layout(annotations=latexannotationshklm)
hklmboot_latex
#Over other inputs
#Capital
hkinpboot <- subplot(hklqboot, hkmqboot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))))
hkinpboot <- hkinpboot %>% layout(annotations=annotationsab)
hkinpboot
#Labor
hlinpboot <- subplot(hlkqboot, hlmqboot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))))
hlinpboot <- hlinpboot %>% layout(annotations=annotationsab)
hlinpboot
#Materials
hminpboot <- subplot(hmkqboot, hmlqboot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
hminpboot <- hminpboot %>% layout(annotations=annotationsab)
hminpboot
#Scale and Ratio Effects
hratio <- subplot(hklboot, hkmboot, hlmboot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-(k-l)", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital/Labor Effect", titlefont=list(size=18), tickfont=list(size=14))),
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-(k-m)", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital/Materials Effect", titlefont=list(size=18), tickfont=list(size=14))),
		scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-(l-m)", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor/Materials Effect", titlefont=list(size=18), tickfont=list(size=14))))
hsboot
hratio
omgboot





