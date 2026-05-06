
hexSticker::sticker(
  #package name
  package="surrosurv",
  # p_family="sans",
  p_size=20,
  p_color="#012b6e",
  p_x=1.05, p_y=0.49,
  #hexagon
  h_fill = "white",
  h_color = "#012b6e",
  h_size = 1.2,
  #subplot
  subplot= "inst/figures/km2.png", #sinon carré qui sort du logo
  s_x=1, s_y=1.15,
  s_width=0.63,
  asp = 1,
  # s_height=5000,
  #output
  filename="inst/figures/logo.png"
)

# magick::image_read('inst/figures/logo.png') |> plot()
