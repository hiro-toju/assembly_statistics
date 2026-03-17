install.packages("devtools")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyr")
install.packages("forcats")  

library(forcats)
library(ggplot2)
library(dplyr)
library(tidyr)

devtools::install_github("hiroakif93/R-functions/AnalysisHelper", force=TRUE)
library(AnalysisHelper)

sml <- readRDS('table/table_ews.rds')


{
  sml |>
  select(everything(), "Observed temporal community change"=`Community change`, "Basin entropy"=`Stable-state entropy`)|>
  pivot_longer(cols=c(`Local structural stability`, `Basin entropy` ))|>
    mutate(
      name = fct_relevel(name, "Local structural stability")
    )|>
  ggplot(aes(x=value, y=`Observed temporal community change`))+
  geom_hline(yintercept=0.5, linetype=4, size=0.3)+
  geom_point(aes(color=as.factor(replicate.id)), size=0.4)+
  facet_wrap(~name, scales="free",   strip.position = "bottom")+
  scale_color_viridis_d()+
    theme_text(unitsize = 0.4)+
    theme(strip.placement = "outside")+
    labs(x="", color="Replicate")
}|>
  ggsave(filename="result/ews.pdf", w=10, h=6, unit="cm")

