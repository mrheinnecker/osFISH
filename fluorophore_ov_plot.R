
test <- read_tsv("/home/rheinnec/Downloads/fluorophores - Sheet1 (1).tsv")





ov_plot_fluorophores <- test %>%
  pivot_longer(
    cols = c("ex", "em"),
    names_to = "cat",
    values_to = "wl"
  ) %>%
  ggplot()+
  geom_segment(
    aes(x=wl, xend=wl, y=0, yend=1, color=fluorophore, linetype=cat)
  )+
  geom_text(
    aes(label = fluorophore,
        x=wl,
        y=1,
        color=info
        ),
    angle=55,hjust=0
  )+
  scale_y_continuous(limits = c(0,1.3))+
  theme_bw()


pdf("/g/schwab/Marco/projects/osFISH/fluorophore_order_plot.pdf", width=15, height=7)
ov_plot_fluorophores
dev.off()






