library('shiny')
source('../R/helpers.R')


ui <- fluidPage(
  titlePanel('Configuration'),
 
  sidebarLayout(
    sidebarPanel(
      selectInput('gtype',
                  'Type of Graph',
                  choices = list(
                    'Erdos-Renyi', 'Watts-Strogatz', 'Power-Law', 'Geometric'
                  )
      ),
      selectInput('network_measure',
                  'Network Measures',
                  choices = list(
                    'Betweenness', 'Closeness', 'Degree', 'Strength', 'Eigenvector'
                  )
      ),
      selectInput('causal_measure',
                  'Causal Measures',
                  choices = list(
                    'KL Causal Effect', 'ACE'
                  )
      ),
      selectInput('estimator',
                  'Estimator',
                  choices = list(
                    'Graphical-Lasso', 'CI-tests (alpha = 0.01)'
                  )
      ),
      sliderInput('nodes', 'Nodes', min = 5, max = 50, value = 10),
      sliderInput('conn', 'Network Density', min = .1, max = .9, value = .5, step = .1),
      sliderInput('edge_weight', 'Average Edge weight', min = 0, max = .5, value = .3, step = .1),
      actionButton('run','Run!')
    ),
 
    mainPanel(
      plotOutput('compplot'),
      plotOutput('corplot')
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  out <- eventReactive(input$run, {
    
    mapping <- list(
      'Watts-Strogatz' = 'watts',
      'Erdos-Renyi' = 'er',
      'Power-law' = 'power',
      'Geometric' = 'geometric'
    )
    
    g <- randDAG(
      n = input$nodes, d = conn2degree(input$nodes, input$conn),
      method = mapping[[input$gtype]], wFUN = list(rnorm, mean = input$edge_weight, sd = .5)
    )
    
    dat <- rDAG(n = 5000, dag = g)
    
    if (input$estimator == 'Graphical-Lasso') {
      G <- get_mrf(dat)
    } else {
      G <- get_mrf2(dat, alpha = .01)
    }
    
    list(g = g, dat = dat, G = G, network_measure = input$network_measure, causal_measure = input$causal_measure)
  })
  
  output$corplot <- renderPlot({
    out <- out()
    network_measure <- ifelse(out$network_measure == 'Degree', 'InDegree', out$network_measure)
    causal_measure <- ifelse(out$causal_measure == 'KL Causal Effect', 'KL', 'ACE')
    plot_single(out$g, out$G, network_measure, causal_measure)
  })
  
  output$compplot <- renderPlot({
    out <- out()
    network_measure <- ifelse(input$network_measure == 'Degree', 'InDegree', input$network_measure)
    causal_measure <- ifelse(out$causal_measure == 'KL Causal Effect', 'KL', 'ACE')
    compplot(out$g, out$G, network_measure, causal_measure)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)