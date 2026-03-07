```@meta
EditURL = "../scripts/lv_ode.jl"
```

````@example lv_literate
using DrWatson
@quickactivate "ProjectoSIR"

using DifferentialEquations
using DataFrames
using StatsPlots
using LaTeXStrings
using Plots
using Statistics
using FFTW

"""
    Modelo de Lotka-Volterra (Presas-Predadores)

    Sistema de equações:
    dx/dt = α*x - β*x*y  # Variação da população de presas
    dy/dt = δ*x*y - γ*y   # Variação da população de predadores

    Onde:
    - x: população de presas (ex: lebres)
    - y: população de predadores (ex: raposas)
    - α: taxa de crescimento natural das presas (na ausência de predadores)
    - β: taxa de predação (eficiência dos predadores em capturar presas)
    - δ: taxa de conversão de presas em novos predadores
    - γ: taxa de mortalidade natural dos predadores (na ausência de presas)
"""
function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, δ, γ = p

    @inbounds begin
        du[1] = α*x - β*x*y  # equação para presas
        du[2] = δ*x*y - γ*y   # equação para predadores
    end
    nothing
end
````

Parâmetros da modelo

````@example lv_literate
p_lv = [
    0.1,   # α: taxa de reprodução das presas
    0.02,  # β: taxa de predação
    0.01,  # δ: eficiência de conversão (presas -> predadores)
    0.3    # γ: mortalidade dos predadores
]
````

Condições iniciais: [presas, predadores]

````@example lv_literate
u0_lv = [40.0, 9.0]
````

Parâmetros de tempo

````@example lv_literate
tspan_lv = (0.0, 200.0)  # tempo total de simulação
dt_lv = 0.01              # passo de integração
````

Criar e resolver o problema

````@example lv_literate
prob_lv = ODEProblem(lotka_volterra!, u0_lv, tspan_lv, p_lv)
sol_lv = solve(prob_lv,
               Tsit5(),              # método de 5ª ordem
               dt = dt_lv,
               reltol = 1e-8,         # tolerância relativa
               abstol = 1e-10,         # tolerância absoluta
               saveat = 0.1,           # salvar a cada 0.1 unidades de tempo
               dense = true)           # permitir interpolação
````

Preparar DataFrame com resultados

````@example lv_literate
df_lv = DataFrame()
df_lv[!, :t] = sol_lv.t
df_lv[!, :prey] = [u[1] for u in sol_lv.u]      # presas
df_lv[!, :predator] = [u[2] for u in sol_lv.u]  # predadores
````

Calcular derivadas para análise

````@example lv_literate
df_lv[!, :dprey_dt] = p_lv[1] .* df_lv.prey .- p_lv[2] .* df_lv.prey .* df_lv.predator
df_lv[!, :dpredator_dt] = p_lv[3] .* df_lv.prey .* df_lv.predator .- p_lv[4] .* df_lv.predator
````

Calcular pontos de equilíbrio

````@example lv_literate
x_star = p_lv[4] / p_lv[3]  # γ/δ - ponto de equilíbrio para presas
y_star = p_lv[1] / p_lv[2]  # α/β - ponto de equilíbrio para predadores

println("="^60)
println("Modelo de Lotka-Volterra (Presas-Predadores)")
println("="^60)
println("\nParâmetros da modelo:")
println("α (taxa de reprodução das presas) = ", p_lv[1])
println("β (taxa de predação) = ", p_lv[2])
println("δ (eficiência de conversão) = ", p_lv[3])
println("γ (mortalidade dos predadores) = ", p_lv[4])

println("\nCondições iniciais:")
println("Presas (x0) = ", u0_lv[1])
println("Predadores (y0) = ", u0_lv[2])

println("\nPontos de equilíbrio:")
println("x* = γ/δ = ", round(x_star, digits=3))
println("y* = α/β = ", round(y_star, digits=3))
````

Gráfico 1: Dinâmica das populações

````@example lv_literate
plt1 = plot(df_lv.t, [df_lv.prey df_lv.predator],
            label = [L"Presas (x)" L"Predadores (y)"],
            xlabel = "Tempo", ylabel = "População",
            title = "Modelo Lotka-Volterra: Dinâmica das Populações",
            linewidth = 2,
            legend = :topright,
            grid = true,
            size = (900, 500),
            color = [:green :red])
````

Adicionar linhas dos pontos de equilíbrio

````@example lv_literate
hline!(plt1, [x_star], color=:green, linestyle=:dash, alpha=0.5, label="x* (equilíbrio presas)")
hline!(plt1, [y_star], color=:red, linestyle=:dash, alpha=0.5, label="y* (equilíbrio predadores)")
````

Salvar gráfico

````@example lv_literate
plots_dir = plotsdir()
if !isdir(plots_dir)
    mkpath(plots_dir)
end
savefig(plt1, joinpath(plots_dir, "lv_dynamics.png"))

println("\nGráfico salvo em: ", joinpath(plots_dir, "lv_dynamics.png"))
println("\nSimulação concluída com sucesso!")
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

