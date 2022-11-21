Base.Base.@kwdef mutable struct state
    x::Float64 = 0.0
    y::Float64 = 1.0
    speedX::Float64 = 1.0
    speedY::Float64 = -0.5
end

function updateWave(timeInterval::Float64, x::Float64, speed::Float64)
    x_next = x + timeInterval*speed

    ## Boundary refelction
    if x_next > 1
        speed *= -1
        x_next = 1 + timeInterval*speed
    elseif x_next < 0
        speed *= -1
        x_next = timeInterval * speed
    end
end



function main()
    state = State()
    const fps::Int64 = 100
    const timeInterval::Float64 = 1.0/fps


    for i in 1:1:100
        state.x = updateWave(state.timeInterval, state.x, state.speedX)
        state.y = updateWave(state.timeInterval, state.y, state.speedY)
    end

end

main()