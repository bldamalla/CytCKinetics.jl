## only for speeding up loading Makie

using CairoMakie

model(x, p) = @. p[1] * x / (x + p[2])

p0 = [60, 0.5]
xs = 0:0.5:5.0; len = length(xs)
xsfine = 0:0.01:5.0
ys = model(xs, p0) .+ randn(len) .* 0.1
ysline = model(xsfine, p0)

f = Figure()

ax = Axis(f[1,1],
    xlabel="Total substrate conc. (μM)",
    ylabel="Initial rate (nM/s)",
    title="Reaction rates vs. concentration"
)

scatter!(ax, xs, ys; label="Data", color=:blue)
lines!(ax, xsfine, ysline; label="Fitted line", color=:red)

yslow = ysline .* 0.95
yshi = ysline .* 1.05
band!(ax, xsfine, yslow, yshi; color=(:red, 0.5))

Legend(f[1,2], ax, framevisible=false)

save("thing.pdf", f)
