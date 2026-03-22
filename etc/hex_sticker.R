require(hexSticker)
require(ggplot2)
require(svglite)

# Create a simple non-linear curve plot
p <- ggplot(data.frame(x = c(-2, 2)), aes(x)) +
  stat_function(fun = function(x) x^3 - x, color = "#ff7f0e", linewidth = 1.5) +
  stat_function(fun = function(x) 1.5*x - 1.7, color = "white", linewidth = 0.8, linetype = "dashed", alpha = 0.7) + # not the actual 1st derivative, but better for the look
  geom_hline(yintercept = 0, linetype = "dashed", color = "white", alpha = 0.7, linewidth = 0.8) +
  geom_point(aes(x = c(1, 1), y = c(0,0)), color = "purple", size = 2) + # Highlighted root
  theme_void() +
  theme_transparent()

# Generate the SVG sticker
s <- hexSticker::sticker(
  p,
  package="nleqslv", p_size=8, p_fontface = "plain", p_family = "mono",
  s_x=1.1, s_y=.75, s_width=1.3, s_height=1,
  h_fill="#1a2b3c", h_color="#cccccc",
  filename="etc/nleqslv_sticker.svg")

plot(s)
