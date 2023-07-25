@setup_workload begin
    truth = [21 9 22 10 23 11 24 12
        17 5 18 6 19 7 20 8
        13 1 14 2 15 3 16 4
        45 33 46 34 47 35 48 36
        41 29 42 30 43 31 44 32
        37 25 38 26 39 27 40 28
        69 57 70 58 71 59 72 60
        65 53 66 54 67 55 68 56
        61 49 62 50 63 51 64 52]
    single = fill(313.0, 35000) .+ randn(35000)
    multi = fill(313.0, (35000, 72)) .+ randn((35000, 72))

    @compile_workload begin
        # Utils
        moving_stats(ones(5000), 3)
        moving_stats(ones(5000), fill(3, 5000))
        mat_shape(ones(72))
        mat_shape(ones(5000, 72))
        sfm(zeros(512))
        active_sfm(randn(35000, 72), 300)

        # Combiners
        combiner(SNR_MAX(), multi)
        combiner(PCC(), multi)
        combiner(PCC2(), multi)
        combiner(EGC(), multi)
        apply2seg(s -> combiner(MRC_PSD(10), s), multi, 300)

        # Occupancy Detection
        occupancy_detection(multi)
        occupancy_detection(multi; max_dist = true)

        # Motion Detection
        move_detect(single; min_samples = 10)
        move_detect(multi; min_samples = 10)
        move_detect(Solei(), single; min_samples = 10)
    end
end
