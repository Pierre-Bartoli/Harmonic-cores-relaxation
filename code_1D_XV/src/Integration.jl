"""Global time of the simulation
"""
const CURRENT_TIME = zeros(Int64, 1)

const init_tabX = zeros(Float64, NPART)
const init_tabV = zeros(Float64, NPART)

const k1X = zeros(Float64, NPART)
const k1V = zeros(Float64, NPART)
const k2X = zeros(Float64, NPART)
const k2V = zeros(Float64, NPART)
const k3X = zeros(Float64, NPART)
const k3V = zeros(Float64, NPART)
const k4X = zeros(Float64, NPART) 
const k4V = zeros(Float64, NPART)

const X1 = zeros(Float64, NPART)
const V1 = zeros(Float64, NPART)
const prev_X1 = zeros(Float64, NPART)
const prev_V1 = zeros(Float64, NPART)

"""
Compute the yoshida algorithm order 1
"""
function YO1!()

    for i=2:NPART
        for j = 1:(i-1)
            pair_two_body_real!(i, j, DT)
        end
    end

    backward_real!(DT)

end

"""
Compute the yoshida algorithm order 2
"""
function YO2!(dt::Float64=DT)

    for i=2:NPART
        for j = 1:(i-1)
            pair_two_body_real!(i, j, dt/2)
        end
    end

    for i = reverse(2:NPART)
        for j = reverse(1:(i-1))
            pair_two_body_real!(i, j, dt/2)
        end
    end

    backward_real!(dt)

end

"""
Compute the yoshida algorithm order 4
"""
function YO4!()

    # magic coef for YO4
    x0 = (-2^(1/3))/(2-2^(1/3))
    x1 = (1.0)/(2.0-2^(1/3))

    YO2!(x1*DT)
    YO2!(x0*DT)
    YO2!(x1*DT)

end

"""
function computing GL2
"""
function GL2!()

    #initial tab filled
    for n = 1:NPART

        init_tabX[n] = tabX[n]
        init_tabV[n] = tabV[n]
    end

    for n = 1:NPART

        #fill W_1 with the assumption that W_(n+1) = W_n
        X1[n] = tabX[n]
        V1[n] = tabV[n]

        #then tabW = W_1 to compute speeds with tab_dot
        tabX[n] = X1[n]
        tabV[n] = V1[n]

    end

    #compute the speeds
    tab_dot!()

    for n = 1:NPART

        #drift
        tabX[n] = drift(init_tabX[n], DT, tabXdot[n])
        tabV[n] = drift(init_tabV[n], DT, tabVdot[n])

        #store the previous W_1
        prev_X1[n] = X1[n]
        prev_V1[n] = V1[n]

        #compute the new W_1
        X1[n] = (tabX[n] + init_tabX[n])/2
        V1[n] = (tabV[n] + init_tabV[n])/2

    end

    #compute the max distance
    d = max_distance(X1, V1, prev_X1, prev_V1)

    i = 0

    #iterate
    while (d > 10e-15 && i <= 50)
        
        i += 1
        
        for n = 1:NPART
            #tabW = W_1 to compute speeds with tab_dot
            tabX[n] = X1[n]
            tabV[n] = V1[n]

        end

        #compute the speeds
        tab_dot!()

        for n = 1:NPART

            #drift
            tabX[n] = drift(init_tabX[n], DT, tabXdot[n])
            tabV[n] = drift(init_tabV[n], DT, tabVdot[n])

            #store prev_W1
            prev_X1[n] = X1[n]
            prev_V1[n] = V1[n]

            #compute new W_1
            X1[n] = (tabX[n] + init_tabX[n])/2
            V1[n] = (tabV[n] + init_tabV[n])/2

        end

        d = max_distance(X1, V1, prev_X1, prev_V1)

    end

end

"""
function that compute the max distance of two points
# arguments
+ X1 cords of the first point
+ prev_X1 cords of the second point
+ V1 cords of the first point
+ prev_V1 cords of the second point

"""
function max_distance(X1::Array{Float64, 1}, V1::Array{Float64, 1}, prev_X1::Array{Float64, 1}, prev_V1::Array{Float64, 1})

    max = 0

    for n = 1:NPART

        if (abs(X1[n] - prev_X1[n]) >= max)

            max = abs(X1[n] -  prev_X1[n])

        end 

        if (abs(V1[n] - prev_V1[n]) >= max)

            max = abs(V1[n] -  prev_V1[n])

        end
    end
    
    return max

end

"""
Function that compute RK1
"""
function RK1!()

    #compute velocities
    tab_dot!()

    #drift
    for n = 1:NPART
        tabX[n] = drift(tabX[n], DT, tabXdot[n])
        tabV[n] = drift(tabV[n], DT, tabVdot[n])
    end

end

"""
Function that compute RK2
"""
function RK2!()
    #compute RK2

    #compute velocities
    tab_dot!()

    #first step half DT
    for n = 1:NPART

        init_tabX[n] = tabX[n]
        init_tabV[n] = tabV[n]

        tabX[n] = drift(tabX[n], DT/2, tabXdot[n])
        tabV[n] = drift(tabV[n], DT/2, tabVdot[n])

    end

    tab_dot!()

    #complete DT
    for n = 1:NPART

        tabX[n] = drift(init_tabX[n], DT, tabXdot[n])
        tabV[n] = drift(init_tabV[n], DT, tabVdot[n])

    end
end

function RK3!()

    for n = 1:NPART

        #store the initial tabW
        init_tabX[n] = tabX[n]
        init_tabV[n] = tabV[n]

    end

    tab_dot!()

    for n = 1:NPART

        #compute k1
        k1X[n] = tabXdot[n]
        k1V[n] = tabVdot[n]

        #move k1 step
        tabX[n] = drift(init_tabX[n], DT/2, k1X[n])
        tabV[n] = drift(init_tabV[n], DT/2, k1V[n])

    end

    tab_dot!()

    for n = 1:NPART

        #compute k2
        k2X[n] = tabXdot[n]
        k2V[n] = tabVdot[n]

        #move k2 step
        tabX[n] = drift(init_tabX[n], DT, -k1X[n] + 2*k2X[n])
        tabV[n] = drift(init_tabV[n], DT, -k1V[n] + 2*k2V[n])

    end

    tab_dot!()

    for n = 1:NPART

        #compute k2
        k3X[n] = tabXdot[n]
        k3V[n] = tabVdot[n]

        #move k2 step
        tabX[n] = drift(init_tabX[n], DT, 1/6*k1X[n] + 2/3*k2X[n] + 1/6*k3X[n])
        tabV[n] = drift(init_tabV[n], DT, 1/6*k1V[n] + 2/3*k2V[n] + 1/6*k3V[n])

    end

end

"""
Function that compute RK4
"""
function RK4!()

    for n = 1:NPART

        #store the initial tabW
        init_tabX[n] = tabX[n]
        init_tabV[n] = tabV[n]

    end

    tab_dot!()

    for n = 1:NPART

        #compute k1
        k1X[n] = tabXdot[n]
        k1V[n] = tabVdot[n]

        #move k1 step
        tabX[n] = drift(init_tabX[n], DT/2, k1X[n])
        tabV[n] = drift(init_tabV[n], DT/2, k1V[n])

    end

    tab_dot!()

    for n = 1:NPART

        #compute k2
        k2X[n] = tabXdot[n]
        k2V[n] = tabVdot[n]

        #move k2 step
        tabX[n] = drift(init_tabX[n], DT/2, k2X[n])
        tabV[n] = drift(init_tabV[n], DT/2, k2V[n])

    end
    
    tab_dot!()

    for n = 1:NPART

        #compute k3
        k3X[n] = tabXdot[n]
        k3V[n] = tabVdot[n]

        #move k3 step
        tabX[n] = drift(init_tabX[n], DT, k3X[n])
        tabV[n] = drift(init_tabV[n], DT, k3V[n])

    end
    
    tab_dot!()

    for n = 1:NPART

        #compute k4
        k4X[n] = tabXdot[n]
        k4V[n] = tabVdot[n]

        #move full step
        tabX[n] = drift(init_tabX[n], DT/6, k1X[n] + 2*k2X[n] + 2*k3X[n] + k4X[n])
        tabV[n] = drift(init_tabV[n], DT/6, k1V[n] + 2*k2V[n] + 2*k3V[n] + k4V[n])

    end

end

if (INTEGRATION == "YO1")

    const integrate_DT!() = YO1!()

elseif (INTEGRATION == "YO2")

    const integrate_DT!() = YO2!()

elseif (INTEGRATION == "YO4")

    const integrate_DT!() = YO4!()

elseif (INTEGRATION == "RK1")

    const integrate_DT!() = RK1!()

elseif (INTEGRATION == "RK2")

    const integrate_DT!() = RK2!()

elseif (INTEGRATION == "RK3")

    const integrate_DT!() = RK3!()

elseif (INTEGRATION == "RK4")

    const integrate_DT!() = RK4!()

elseif (INTEGRATION == "GL2")

    const integrate_DT!() = GL2!()

else

    const integrate_DT!() = RK2!()

end

"""
function that make a drift
# arguments
+ init_tab, which the start point
+ dt
+ drift_value, is the value that is made by the drift_value
# return 
+ the drifted value

"""
function drift(init_tab, dt, drift_value)

    return init_tab + dt*drift_value

end
"""
Wrapped function that performs the integration
for NSTEPS steps
"""
function integrate_NSTEPS!()

    for istep=1:(NSTEPS) # Loop over the steps to perform

        integrate_DT!() # Integrating for one timestep
        CURRENT_TIME[1] += 1 #update the time
    end

end

"""
Give the current time of the simulation
# Return 
+ current time

"""
function get_time()
    return CURRENT_TIME[1]*DT
end