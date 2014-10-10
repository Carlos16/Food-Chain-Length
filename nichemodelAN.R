web.files <- list.files(path = "C:/Users/jjborrelli/Desktop/CHAINDATA_ni/", pattern = "web")
temp.ls <- list()
for(i in 1:length(web.files)){
  temp.ls[[i]]<- fread(paste("C:/Users/jjborrelli/Desktop/CHAINDATA_ni/", 
                             web.files[i], sep = ""))
}
webdata <- rbindlist(temp.ls)
webdata$C <- factor(webdata$C)

sem.l <- function(x){mean(x) - 1.96*sqrt(var(x)/length(x))}
sem.u <- function(x){mean(x) + 1.96*sqrt(var(x)/length(x))}

sub1 <- subset(webdata, par == "1/-10" | par == "10/-1" | par == "1/-5" | par == "5/-1")
sub2 <- subset(webdata, par == "5/-10" | par == "10/-5")
sub3 <- subset(webdata, par == "1/-1" | par == "10/-10" | par == "5/-5")

ggplot(sub1, aes(x = factor(diam+1), y = qss)) + 
  geom_point(alpha = .25, position = position_jitter(w=0.2), col = "grey58") + 
  stat_summary(fun.y="mean", geom="point") +
  stat_summary(fun.ymin = sem.l, fun.y = "mean", fun.ymax = sem.u,
               geom="errorbar", width = .2) + 
  #geom_smooth(method = "glm", aes(lty = C)) +
  facet_grid(par~C) + theme_bw() +
  xlab("Longest Food Chain Length") + ylab("Quasi sign-stability")

ggsave("C:/Users/jjborrelli/Desktop/CHAINDATA_ni/niche1.jpeg", dpi = 1000, width = 8, height = 10)
dev.off()

ggplot(sub2, aes(x = factor(diam+1), y = qss)) + 
  geom_point(alpha = .25, position = position_jitter(w=0.2), col = "grey58") + 
  stat_summary(fun.y="mean", geom="point") +
  stat_summary(fun.ymin = sem.l, fun.y = "mean", fun.ymax = sem.u,
               geom="errorbar", width = .2) + 
  #geom_smooth(method = "glm", aes(lty = C)) +
  facet_grid(par~C) + theme_bw() +
  xlab("Longest Food Chain Length") + ylab("Quasi sign-stability")

ggsave("C:/Users/jjborrelli/Desktop/CHAINDATA_ni/niche2.jpeg", dpi = 1000, width = 8, height = 10)
dev.off()

ggplot(sub3, aes(x = factor(diam+1), y = qss)) + 
  geom_point(alpha = .25, position = position_jitter(w=0.2), col = "grey58") + 
  stat_summary(fun.y="mean", geom="point") +
  stat_summary(fun.ymin = sem.l, fun.y = "mean", fun.ymax = sem.u,
               geom="errorbar", width = .2) + 
  #geom_smooth(method = "glm", aes(lty = C)) +
  facet_grid(par~C) + theme_bw() +
  xlab("Longest Food Chain Length") + ylab("Quasi sign-stability")

ggsave("C:/Users/jjborrelli/Desktop/CHAINDATA_ni/niche3.jpeg", dpi = 1000, width = 8, height = 10)
dev.off()

ggplot(webdata, aes(x = factor(diam+1), y = qss)) + 
  geom_point(alpha = .25, position = position_jitter(w=0.2), col = "grey58") + 
  stat_summary(fun.y="mean", geom="point") +
  stat_summary(fun.ymin = sem.l, fun.y = "mean", fun.ymax = sem.u,
               geom="errorbar", width = .2) + 
  #geom_smooth(method = "glm", aes(lty = C)) +
  facet_grid(par~C) + theme_bw() +
  xlab("Longest Food Chain Length") + ylab("Quasi sign-stability")

ggsave("C:/Users/jjborrelli/Desktop/CHAINDATA_ni/niche-all.jpeg", dpi = 1000, width = 8, height = 10)
dev.off()


df <- data.frame(C = sub1$C, D = factor(sub1$diam +1), qss = sub1$qss, par = factor(sub1$par))
df2 <- data.frame(C = sub2$C, D = factor(sub2$diam +1), qss = sub2$qss, par = factor(sub2$par))
df3 <- data.frame(C = sub3$C, D = factor(sub3$diam +1), qss = sub3$qss, par = factor(sub3$par))


ggplot(aes(y = qss, x = D), data = df) + geom_boxplot() + 
  facet_grid(par~C, scale = "fixed") + theme_bw() +
  xlab("Longest Chain Length") + ylab("Quasi Sign-Stability")

ggplot(aes(y = qss, x = D), data = df2) + geom_boxplot() +
  facet_grid(par~C, scale = "fixed") + theme_bw() + 
  xlab("Longest Chain Length") + ylab("Quasi Sign-Stability")

ggplot(aes(y = qss, x = D), data = df3) + geom_boxplot() + 
  facet_grid(par~C, scale = "fixed") + theme_bw() + 
  xlab("Longest Chain Length") + ylab("Quasi Sign-Stability")