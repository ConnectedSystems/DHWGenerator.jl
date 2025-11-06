using Random
using OrderedCollections

import GeoDataFrames as GDF
using DataFrames

using Distances
import GeometryOps as GO

using YAXArrays
using NetCDF
using NCDatasets

"""
    generate_dhw_trajectories(
        n_years::Int64,
        start_year::Int64,
        end_year::Int64;
        rng::AbstractRNG=Random.GLOBAL_RNG,
        with_dhw=true,
        warming_rate::Float32=0.15f0,
        seasonal_amplitude::Float32=1.2f0,
        dhw_threshold::Float32=4.0f0,
        spatial_length_scale::Float32=0.5f0,
        noise_amplitude::Float32=0.9f0
    )

Generate synthetic environmental data for coral reef modeling, specifically Degree Heating
Weeks (DHW) trajectories that simulate realistic marine heatwave patterns under climate
change scenarios.

*TODO: A key limitation currently is that there is very limited consideration of the
spatial correlation. For that, we would need some measure of how close sites/reefs are to
each other.*

# Extended help

Generate plausible DHW time series by combining multiple environmental components:

1. **Long-term warming trend**: Simulates gradual ocean warming (like RCP4.5 scenarios)
2. **Seasonal cycles**: Models natural temperature variations throughout the year
3. **Weather noise**: Adds realistic short-term temperature fluctuations
4. **Spatial variation**: Creates location-specific temperature offsets
5. **Acute heatwave events**: Superimposes extreme marine heatwave events
6. **DHW accumulation**: Converts temperature anomalies to ecologically-relevant DHW values

The DHW calculation follows coral bleaching research where:
- DHW accumulates when temperatures exceed a threshold (default 4°C above baseline)
- Values decay over time when temperatures drop
- Extreme events can cause rapid DHW spikes that lead to mass bleaching

## Temperature Anomaly Construction
For each location and timestep, the temperature anomaly is built from:

```julia
temp_anomaly = warming_trend + seasonal_cycle + weather_noise + spatial_offset
```

Where:
- **warming_trend**: Linear increase over time (`warming_rate * years_elapsed`)
- **seasonal_cycle**: Sinusoidal pattern with amplitude that increases over time
- **weather_noise**: Random normal variations that get larger in later years
- **spatial_offset**: Location-specific random offset (0 - 0.8°C)

## DHW Accumulation Rules
- **Above threshold**: DHW accumulates as `(temp_anomaly - threshold) / 4.0`
- **Below threshold**: DHW decays rapidly (`previous_DHW * 0.7`)
- **During heating**: DHW decays slowly (`previous_DHW * 0.92`)
- **Soft cap**: Values above 20 DHW are dampened but can still fluctuate

## Acute Event Generation
The function adds realistic extreme heatwave events:
- **Frequency**: ~1 event per year on average (`n_timesteps / 12`)
- **Probability**: Increases over time (simulating worsening climate)
- **Duration**: 2-5+ weeks, longer in later years
- **Intensity**: 8-25+ DHW, stronger in later years
- **Shape**: Rapid onset (30% of duration) → peak → gradual decline (30% of duration)

## Ecological Realism
The parameters are tuned based on coral bleaching research:
- **4 DHW**: Threshold where bleaching typically begins
- **8+ DHW**: Significant bleaching and mortality expected
- **20+ DHW**: Severe bleaching events (soft-capped with fluctuations)

# Arguments
- `n_years`: Number of time steps to simulate (in years)
- `start_year`: Starting year for the simulation
- `end_year`: End year for the simulation
- `rng`: Random number generator for reproducible results
- `warming_rate`: Rate of long-term warming per year in °C (default: 0.15°C/year)
- `seasonal_amplitude`: Strength of seasonal temperature cycles (default: 1.2°C)
- `dhw_threshold`: Temperature threshold above which DHW accumulates (default: 4.0°C)
- `spatial_length_scale`: Controls how strongly nearby reef sites influence each other - adjust as needed (default: 0.5)
- `noise_amplitude`: Magnitude of random weather variations (default: 0.9°C)
- `avg_extreme_years`: Average time between extreme years (default: 12 years)
"""
function generate_dhw_trajectories(
    reef_representation::String,
    start_year::Int64,
    end_year::Int64;
    rng::AbstractRNG=Random.GLOBAL_RNG,
    warming_rate::Float32=0.15f0,
    seasonal_amplitude::Float32=1.2f0,
    dhw_threshold::Float32=4.0f0,
    spatial_length_scale::Float32=0.5f0,
    noise_amplitude::Float32=0.9f0,
    avg_extreme_years::Int64=12
)::YAXArray
    n_years = (end_year - start_year) + 1

    # Read in geospatial data
    geo_coord = GDF.read(reef_representation)
    n_locations = nrow(geo_coord)

    # Generate centroids of spatial polygons and extract latitudes and longitudes
    cents = GO.centroid.(geo_coord.geometry)
    longs = first.(cents)
    lats = last.(cents)
    coords = hcat(longs, lats)

    # Initialize vectors
    dhw_data = zeros(Float32, n_years, n_locations)
    spatial_kernel = zeros(Float32, n_locations, n_locations)

    for i in 1:n_locations, j in 1:n_locations
        # Compute pairwise haversine distances (in kilometers)
        c = haversine(coords[i, :], coords[j, :]) * 0.001

        # Compute spatial correlation using Gaussian spatial kernel
        spatial_kernel[i, j] = Float32(exp(-c^2 / (2.0 * spatial_length_scale^2)))
    end

    # Calculate baseline parameters
    years = range(start_year, length=n_years) .- start_year

    for loc in 1:n_locations
        # Location-specific random offset to create spatial variation
        spatial_offset = rand(rng, Float32) * 0.8f0  # Increased spatial variation

        # Track consecutive weeks above threshold
        weeks_above_threshold = 0

        # Generate temperature anomaly time series
        for t in 1:n_years
            # Generate spatially correlated noise for this timestep
            base_noise = randn(rng, n_locations)
            spatially_correlated_noise = spatial_kernel * base_noise

            # Long-term warming trend
            warming_trend = warming_rate * years[t]

            # Seasonal cycle with increasing amplitude over time
            seasonal_scaling = 1.0f0 + 0.2f0 * (years[t] / years[end])
            seasonal_cycle = seasonal_amplitude * seasonal_scaling * sin(2π * (t % 12) / 12)

            # Random weather variations with increased variability over time
            weather_scaling = 1.0f0 + 0.3f0 * (years[t] / years[end])
            weather_noise = noise_amplitude * weather_scaling * spatially_correlated_noise[loc]

            # Combined temperature anomaly
            temp_anomaly = (
                warming_trend +
                seasonal_cycle +
                weather_noise +
                spatial_offset
            )

            prev_dhw = t > 1 ? dhw_data[t-1, loc] : 0.0f0

            # Convert temperature anomaly to DHW with potential acute events
            if temp_anomaly > dhw_threshold
                weeks_above_threshold += 1

                # Base DHW accumulation (where 4.0 converts weekly temperature to DHW)
                # Assumes 1 DHW = 1°C > threshold for 1 week
                dhw_accumulation = (temp_anomaly - dhw_threshold) / 4.0f0

                # Add possibility of acute temperature spikes
                if weeks_above_threshold >= 2
                    # Chance of acute event increases with warming trend
                    acute_probability = min(0.2f0 * (1.0f0 + warming_trend), 0.4f0)

                    if rand(rng, Float32) < acute_probability
                        # Generate acute spike with magnitude increasing over time
                        time_factor = years[t] / years[end]
                        base_spike = 3.0f0 + 2.0f0 * time_factor  # Spikes get larger over time
                        spike_magnitude = rand(rng) * base_spike + 2.0f0
                        dhw_accumulation += spike_magnitude
                    end
                end

                # Calculate new DHW with modified decay
                # Reduce value for slower decay during heating
                base_dhw = prev_dhw * 0.92f0
                raw_dhw = base_dhw + dhw_accumulation

                # Apply soft cap with fluctuations
                soft_cap = 20.0f0
                if raw_dhw > soft_cap
                    excess = raw_dhw - soft_cap
                    damping_factor = 1.0f0 / (1.0f0 + 0.3f0 * excess)

                    # Larger fluctuations in later years
                    time_factor = years[t] / years[end]
                    max_fluctuation = 4.0f0 + 2.0f0 * time_factor
                    fluctuation = (rand(rng) - 0.5f0) * min(max_fluctuation, excess * 0.6f0)

                    dhw_data[t, loc] = soft_cap + (excess * damping_factor) + fluctuation
                else
                    dhw_data[t, loc] = raw_dhw
                end
            else
                # Reset counter and apply faster decay when below threshold
                weeks_above_threshold = 0
                dhw_data[t, loc] = max(0.0f0, prev_dhw * 0.7f0)
            end
        end

        # Add extreme marine heatwave events
        n_extreme_events = floor(Int, n_years / avg_extreme_years)

        for _ in 1:n_extreme_events
            event_time::Int64 = rand(rng, 1:n_years)
            time_progress = years[event_time] / years[end]
            event_probability = time_progress * 1.8f0  # Higher probability in later years

            if rand(rng, Float32) < event_probability
                # Duration increases with time
                base_duration = 2:5
                extra_duration = rand(rng, 0:floor(Int, 3 * time_progress))
                event_duration = rand(rng, base_duration) + extra_duration

                # Base magnitude increases with time (8 - 17 DHW)
                base_magnitude = 8.0f0 + 17.0f0 * time_progress
                event_magnitude = rand(rng) * 5.0f0 + base_magnitude

                # Add the extreme event with realistic onset/decline
                for t in event_time:min(event_time + event_duration, n_years)
                    relative_pos = (t - event_time) / event_duration
                    if relative_pos <= 0.3f0
                        # Rapid onset
                        scaling = relative_pos / 0.3f0
                    elseif relative_pos >= 0.7f0
                        # Gradual decline
                        scaling = 1.0f0 - ((relative_pos - 0.7f0) / 0.3f0)
                    else
                        # Peak
                        scaling = 1.0f0
                    end

                    dhw_value = event_magnitude * scaling
                    dhw_data[t, loc] = max(dhw_data[t, loc], dhw_value)
                end
            end
        end
    end

    # Variable axes
    v_axes = (
        Dim{:timestep}(start_year:end_year),
        Dim{:location}(geo_coord.UNIQUE_ID)
    )
    prop = OrderedDict(
        "longitude" => longs,
        "latitude" => lats
    )
    dhw = YAXArray(v_axes, dhw_data, prop)

    return dhw
end

function write_dataset(fn::String, dhw_dataset::YAXArray)::Nothing

    try
        # Create NetCDF file
        ds = NCDatasets.Dataset(fn, "c")

        # Define dimensions
        defDim(ds, "time", size(dhw_dataset, 1))
        defDim(ds, "location", size(dhw_dataset, 2))

        # Define variables and attributes
        lat_var = defVar(ds, "latitude", Float32, ("location",), attrib=OrderedDict(
            "units" => "degrees_south",
            "standard_name" => "latitude",
            "long_name" => "latitude",))
        lon_var = defVar(ds, "longitude", Float32, ("location",), attrib=OrderedDict(
            "units" => "degrees_east",
            "standard_name" => "longitude",
            "long_name" => "longitude",))
        dhw_var = defVar(ds, "DHW", Float32, ("time", "location"), attrib=OrderedDict(
            "units" => "DegC-weeks",
            "standard_name" => "DHW",
            "long_name" => "degree heating weeks",))

        #Write data
        lat_var[:] = dhw_dataset.properties["latitude"]
        lon_var[:] = dhw_dataset.properties["longitude"]
        dhw_var[:, :] = dhw_dataset.data

        close(ds)
    catch err
        @info "Error encountered when attempting to create netCDF: $(err)"
        close(ds)
    end

    return nothing
end
