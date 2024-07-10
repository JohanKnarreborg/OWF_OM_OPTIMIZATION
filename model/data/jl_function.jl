


function wind_to_power(wind_speeds, rated_power, rated_wind_speed, cut_in_speed, cut_out_speed, area, power_coefficient)
    # Constants
    air_density = 1.225  # kg/m^3, standard density of air at sea level
    # Initialize an array to store power generated for each wind speed
    power_generated = Float64[]
    # Calculate power for each wind speed
    for v in wind_speeds
        if v < cut_in_speed || v > cut_out_speed
            push!(power_generated, 0.0)
        elseif v <= rated_wind_speed
            # Power proportional to the cube of the wind speed within operational range
            power = 0.5 * air_density * area * (v ^ 3) * power_coefficient
            push!(power_generated, min(power, rated_power))
        else
            # Constant power output at rated power once rated wind speed is exceeded
            push!(power_generated, rated_power)
        end
    end
    return power_generated
end