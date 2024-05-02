library(shiny)
library(ggplot2)
library(gridExtra) # For arranging plots
library(muhaz)
library(survHE)
library(survminer)
library(landest)

load('./shinyAppData.Rdata')

# Define UI
ui <- fluidPage(
  titlePanel("Blended Model Visualization"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("t1", "Start of blending interval:", min = 14, max = 24, value = 24),
      sliderInput("t2", "End of blending interval:", min = 36, max = 60, value = 60),
      sliderInput("a", "Beta parameter a:", min = 0.1, max = 10, value = 5),
      sliderInput("b", "Beta parameter b:", min = 0.1, max = 10, value = 5)
    ),
    mainPanel(
      plotOutput("hazPlot"),
      plotOutput("survPlot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  output$hazPlot <- renderPlot({
    blended_model(input$t1, input$t2, input$a, input$b, plot_type = "haz")
  })
  
  output$survPlot <- renderPlot({
    blended_model(input$t1, input$t2, input$a, input$b, plot_type = "surv")
  })
}

# Define function to return plots
blended_model <- function(t1, t2, a, b, plot_type) {
  ## Set blending interval and Beta parameters for each arm
  blending_interval_Pem <- list(min = t1, max = t2)
  beta_params_Pem <- list(alpha = a, beta = b)
  
  blending_interval_Ipi <- list(min = t1, max = t2)
  beta_params_Ipi <- list(alpha = a, beta = b)
  
  wt_par_Pem <- list(a=blending_interval_Pem$min, b=blending_interval_Pem$max,
                     shape1=beta_params_Pem$alpha, shape2=beta_params_Pem$beta)
  wt_par_Ipi <- list(a=blending_interval_Ipi$min, b=blending_interval_Ipi$max,
                     shape1=beta_params_Ipi$alpha, shape2=beta_params_Ipi$beta)
  
  ## Weight function for each arm
  
  tp <- seq(0,th,0.05)
  weight_Pem <- with(wt_par_Pem,
                     pbeta((tp - a)/(b - a), shape1, shape2))
  weight_Ipi <- with(wt_par_Ipi,
                     pbeta((tp - a)/(b - a), shape1, shape2))
  
  # Blend hazards
  h_blendh_Pem <- (1 - weight_Pem) * h_chosen_OS_Pem$est + weight_Pem * h_chosen_OS_Scha$est
  h_blendh_Ipi <- (1 - weight_Ipi) * h_chosen_OS_Ipi$est + weight_Ipi * h_chosen_OS_Scha$est
  h_blendh_Ipi[1] <- 0
  
  # Create plots
  if (plot_type == "haz") {
    # Plot blended hazard for Pembrolizumab
    p1 <- ggplot() +
      geom_line(aes(haz_OS_Pem_7y$est.grid, haz_OS_Pem_7y$haz.est, colour="Updated data"), linewidth=1, linetype='dashed') +
      geom_line(aes(tp, h_blendh_Pem, colour="Blended hazard"), linewidth = 1) +
      geom_line(data=subset(h_chosen_OS_Scha, time>mf), aes(x=time, y=est, colour="TA366 base case"), linewidth=1, linetype='dashed') +
      geom_vline(xintercept=blending_interval_Pem$min, linetype="dashed", color="grey") +
      geom_vline(xintercept=blending_interval_Pem$max, linetype="dashed", color="grey") +
      xlab('Time (months)') + ylab('Hazard') +
      scale_x_continuous(limits=c(0,th), breaks=seq(0,th,12)) +
      scale_colour_manual(name="", values=c("Blended hazard"="#F8766D", "Updated data"="black", "TA366 base case"="#619CFF")) +
      theme(legend.direction="vertical", legend.box="vertical", legend.position=c(0.8,0.8)) +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), text=element_text(size=12)) +
      labs(title=paste('Blending interval (', blending_interval_Pem$min, ',', blending_interval_Pem$max, ')', ', Beta(', beta_params_Ipi$alpha, ',', beta_params_Ipi$beta,')', sep=""), subtitle="Pembrolizumab") 
    
    # Plot blended hazard for Ipilumumab
    p2 <- ggplot() +
      geom_line(aes(haz_OS_Ipi_7y$est.grid, haz_OS_Ipi_7y$haz.est, colour="Updated data"), linewidth=1, linetype='dashed') +
      geom_line(aes(tp, h_blendh_Ipi, colour="Blended hazard"), linewidth=1) +
      geom_line(data=subset(h_chosen_OS_Scha,time>mf), aes(x=time, y=est, colour="TA366 base case"),
                linewidth=1, linetype="dashed") +
      geom_vline(xintercept=blending_interval_Ipi$min, linetype="dashed", color="grey") +
      geom_vline(xintercept=blending_interval_Ipi$max, linetype="dashed", color="grey") +
      xlab('Time (months)') + ylab('Hazard') +
      scale_x_continuous(limits=c(0,th), breaks=seq(0,th,12)) +
      scale_colour_manual(name="", values=c("Blended hazard"="#F8766D", "Updated data"="black", "TA366 base case"="#619CFF")) +
      theme(legend.direction="vertical", legend.box="vertical", legend.position=c(0.8,0.8)) +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), text=element_text(size=12)) +
      labs(title=paste('Blending interval (', blending_interval_Ipi$min, ',', blending_interval_Ipi$max, ')', ', Beta(', beta_params_Ipi$alpha, ',', beta_params_Ipi$beta,')', sep=""), subtitle="Ipilumumab")
    
    haz_plot <- ggarrange(p1, p2, nrow=1)
    
    return(haz_plot)
  } else if (plot_type == "surv") {
    # Blended survival 
    S_blendh_Pem <- exp(-cumsum(h_blendh_Pem) * 0.05)
    S_blendh_Ipi <- exp(-cumsum(h_blendh_Ipi) * 0.05)
    
    # Plot survival curves
    surv_plot <- ggplot() + 
      geom_line(aes(x=km_OS_7y_Pem$time, y=km_OS_7y_Pem$surv, linetype="Pembrolizumab", color='Updated Kaplan-Meier'), linewidth=1) +
      geom_line(aes(x=km_OS_7y_Ipi$time, y=km_OS_7y_Ipi$surv, linetype="Ipilimumab", color='Updated Kaplan-Meier'), linewidth=1) +
      geom_line(data=data.frame(tp,S_blendh_Pem), aes(tp, S_blendh_Pem, color="Blended survival", linetype="Pembrolizumab"), linewidth=1) +
      geom_line(data=data.frame(tp,S_blendh_Ipi), aes(tp, S_blendh_Ipi, color="Blended survival", linetype="Ipilimumab"), linewidth=1) +
      geom_line(data=data.frame(t_comp, S_comp_Pem), aes(t_comp,S_comp_Pem, color='TA366 base case', linetype="Pembrolizumab"), linewidth=1) +
      geom_line(data=data.frame(t_comp, S_comp_Pem), aes(t_comp, S_comp_Ipi, color='TA366 base case', linetype="Ipilimumab"), linewidth=1) +
      geom_vline(xintercept=blending_interval_Pem$min, linetype="dashed", color="grey") +
      geom_vline(xintercept=blending_interval_Pem$max, linetype="dashed", color="grey") +
      scale_x_continuous(name="Time (months)", breaks=seq(0,84,12), limits=c(0,84)) +
      scale_y_continuous(name="Overall survival", limits=c(0,1), labels=scales::percent) +
      scale_colour_manual(name="Model", values=c("Blended survival"="#F8766D", "Updated Kaplan-Meier"="black", "TA366 base case"="#619CFF")) +
      scale_linetype_manual(name="Treatment arm", values=c("Pembrolizumab"="solid", "Ipilimumab"="dashed")) +
      theme_classic() +
      theme(legend.direction="vertical", legend.box="vertical", legend.position=c(0.8,0.8)) +
      labs(title=paste('Blending interval (', blending_interval_Pem$min, ',', blending_interval_Pem$max, ')', ', Beta(', beta_params_Ipi$alpha, ',', beta_params_Ipi$beta,')', sep=""))
    
    return(surv_plot)
  }
}

# Run the application 
shinyApp(ui = ui, server = server)
