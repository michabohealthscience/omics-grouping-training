#probability density function R - t-statistic

x <- seq(-5, 5, length.out = 100)
df <- 10
pdf <- dt(x, df)
plot(x, pdf, type="l", main="PDF of t-distribution with 10 degrees of freedom")