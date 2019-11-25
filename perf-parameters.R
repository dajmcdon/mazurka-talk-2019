#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dpf)
library(knitr)
library(tidyverse)
library(RColorBrewer)
library(viridis)
green = '#00884E'
blue = '#053B64'
red = '#FF4900'
orange = '#FF9200'
db = "#053B64"
load("../../Manuscripts/Chopin Mazurka/dpf/extras/mazurkaResults.Rdata") 
source("../../Manuscripts/Chopin Mazurka/dpf/manuscript/dirichlet_precision.R")

pvec_ml = pvec_ml %>% select(-value,-fevals,-gevals,-convergence) %>%
    data.matrix %>% data.frame
th = theme_minimal(base_size = 24,base_family = "Times") +
    theme(axis.text=element_text(color = db),
          legend.title = element_blank(),
          legend.position = 'bottom',
          plot.title=element_text(color = db),
          strip.text = element_text(hjust=0,color=db),
          text = element_text(color = db))
#theme_set(theme_minimal(base_size=15,base_family="Times"))
data(tempos)
lt = diff(c(tempos$note_onset,61))
performers = gsub('_',' ',rownames(pvec_ml))
rownames(pvec_ml) = performers
temp = tempos
colnames(temp)[-c(1:3)] = performers
perfshapes = 15:18
perfcols = viridis_pal(begin=.2)(4)
tt = temp %>% select(-meas_num, -beat) %>% 
    pivot_longer(-note_onset)
hc_parm = pvec_ml %>% Dmats %>% Reduce(f='+')
rownames(hc_parm) <- colnames(hc_parm) <- performers
othercut = .35
subs = apply(hc_parm,1,quantile,probs=4/46) < othercut
sDmat = hc_parm[subs,subs]
nclusts = 4
sdends = sDmat %>% as.dist %>% hclust %>% as.dendrogram
clustered = data.frame(clust = as.factor(cutree(as.hclust(sdends), k = nclusts)),
                       performer = row.names(sDmat))
tt = full_join(tt, clustered,by=c("name"="performer")) 
levels(tt$clust) = c(levels(tt$clust),"other")
tt$clust[is.na(tt$clust)] = "other"

hc_parm = hc_parm %>% as_tibble(rownames = "performer") %>%
    pivot_longer(-performer)
subs = as.vector(lower.tri(matrix(NA,length(performers),length(performers)),diag = TRUE))
hc_parm = hc_parm[subs,]





# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("perf",
                        "Performance",
                        performers, selected="Richter 1976"),
            selectInput("param",
                        "Estimated parameters", performers,
                        selected="Richter 1976")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    ),
    splitLayout(
        plotOutput("perfPlot"),
        plotOutput("parmPlot")
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    plotdat <- reactive({
        y = matrix(unlist(select(temp, input$perf)), nrow=1)
        params = unlist(pvec_ml[input$param,])
        pmats = musicModel(lt, params[1], params[2:4], c(params[5],1,1),
                       params[6:12], c(132,0), c(400,10))
        beam = with(pmats, beamSearch(a0, P0, c(1,0,0,0,0,0,0,0,0,0,0), dt, ct, Tt, Zt,
                                  HHt, GGt, y, transMat, 400))
        bestpath = beam$paths[which.max(beam$weights),]
        kal = kalman(pmats, bestpath, y)
        plotdat = data.frame(measure = tempos$note_onset, tempo = c(y), 
                         inferred = c(kal$ests), state = convert10to4(bestpath))
    
        plotdat$performer = input$perf
        
        plotdat %>% mutate(statesused = state, 
                           state=factor(state,levels=1:4,labels=c('constant','decel','accel','stress')))
    })
    
    output$statePlot <- renderPlot({
        vec = sort(unique(plotdat()$statesused))
        ggplot(plotdat()) + 
            geom_rect(data=data.frame(xmin = 33, xmax = 45, ymin = -Inf, ymax = Inf),
                      aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                      fill = 'gray90', color = 'gray90') +
            geom_line(aes(x=measure, y=tempo), color='gray40') +
            geom_point(aes(x=measure, y=inferred, color=state,shape=state),size=2.5) +
            # scale_color_brewer(palette='Set1') +
            scale_color_manual(values=perfcols[vec]) +
            scale_shape_manual(values=perfshapes[vec]) +
            th +
            facet_wrap(~performer)
    })
    
    output$perfPlot <- renderPlot({
        ggplot(tt) + 
            geom_rect(data=data.frame(xmin = 33, xmax = 45, ymin = -Inf, ymax = Inf),
                      aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                      fill = 'gray90', color = 'gray90') +
            geom_line(aes(y=value,x=note_onset,color=name)) + th +
            scale_color_manual(values = rep(fivecolors,length.out=nrow(pvec_ml))) +
            ylab("tempo (bpm)") + xlab("measure") + 
            theme(legend.position = "none")
    })
    
    output$parmPlot <- renderPlot({
        vv = filter(hc_parm, performer==input$perf, name==input$param)
        ggplot(hc_parm, aes(x=value)) + geom_density(fill=green,alpha=.6) +
            geom_rug() + th + xlab('parameter distance') +
            geom_vline(xintercept = vv$value, color=db, size=2)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
