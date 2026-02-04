using Test
using Catlab
using Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.Programs

# ------------------------------------------------------------
# Julia 1.11-safe stdout capture: redirect_stdout accepts Pipe, not IOBuffer.
# ------------------------------------------------------------
"""
    capture_stdout(f) -> String

Run `f()` while capturing everything printed to STDOUT, return captured text.
Works on Julia 1.11 (uses Pipe()).
"""
function capture_stdout(f::Function)::String
    p = Pipe()
    buf = IOBuffer()

    reader = @async begin
        while !eof(p)
            write(buf, readavailable(p))
            yield()
        end
    end

    redirect_stdout(p) do
        f()
    end

    closewrite(p)
    wait(reader)
    return String(take!(buf))
end


@testset "id_cf_step41!eg1" begin
    @present CausalR_eg1(FreeSymmetricMonoidalCategory) begin
        (R, X, W, Z, Y, D)::Ob
        cUR::Hom(munit(), R)
        cUX::Hom(R, X)
        cUW::Hom(X, W)
        cUY::Hom(R ⊗ W ⊗ Z, Y)
        cUD::Hom(munit(), D)
        cUZ::Hom(D, Z)
        doX::Hom(munit(), X)
        doD_d::Hom(munit(), D)
        doZ::Hom(munit(), Z)
        obsX_x::Hom(X, munit())
        obsD_d::Hom(D, munit())
        obsZ_z::Hom(Z, munit())
    end

    prog_overall = @program CausalR_eg1 () begin
        r  = cUR()
        x_obs = cUX(r)
        obsX_x(x_obs)

        x_do  = doX()
        w     = cUW(x_do)

        z_do  = doZ()
        y     = cUY(r, w, z_do)

        return y
    end

    display_var = Set([:Z, :X, :W, :Y])

    expected_lines = [
        "---- Step 4.1 on R-fragment #1 ----",
        "Fragment 1, var Z, case2 = case2b_rewritten",
        "Fragment 1, var W, case2 = case2a_do_nothing",
        "Fragment 1, var X, case1 = case1b_detected",
        "Fragment 1, var Y, case1 = case1a_do_nothing",
        "---- Step 4.1 on R-fragment #2 ----",
        "Fragment 2, var W, case1 = case1a_do_nothing",
        "Fragment 2, var X, case2 = case2b_rewritten",
    ]

    out = capture_stdout() do
        id_cf_step41!(prog_overall, display_var; verbose=true)
    end

    for line in expected_lines
        @test occursin(line, out)
    end
end


@testset "id_cf_step41!eg2" begin
    @present CausalR_eg2(FreeSymmetricMonoidalCategory) begin
        (R, X, W, Z, Y, D)::Ob
        cUR::Hom(munit(), R)
        cUX::Hom(R ⊗ Y, X)
        cUW::Hom(X, W)
        cUY::Hom(R ⊗ W ⊗ Z, Y)
        cUD::Hom(munit(), D)
        cUZ::Hom(D, Z)
        doX::Hom(munit(), X)
        doD_d::Hom(munit(), D)
        doZ::Hom(munit(), Z)
        doY::Hom(munit(), Y)
        obsX_x::Hom(X, munit())
        obsD_d::Hom(D, munit())
        obsZ_z::Hom(Z, munit())
    end

    prog_overall = @program CausalR_eg2 () begin
        r  = cUR()
        y0 = doY()
        x_obs = cUX(r, y0)
        obsX_x(x_obs)

        x_do  = doX()
        w     = cUW(x_do)

        z_do  = doZ()
        y1    = cUY(r, w, z_do)

        return y1
    end

    display_var = Set([:Z, :X, :W, :Y])

    target = "FAIL (case1, var=Y): does not satisfy 1.a or 1.b conditions"
    err_ref = Ref{Any}(nothing)

    out = capture_stdout() do
        try
            id_cf_step41!(prog_overall, display_var; verbose=true)
        catch e
            err_ref[] = e
        end
    end

    @test err_ref[] !== nothing
    msg = sprint(showerror, err_ref[])
    @test occursin(target, out) || occursin(target, msg)
end


@testset "id_cf_step41!eg3" begin
    @present CausalR_eg3(FreeSymmetricMonoidalCategory) begin
        (R, X, W, Z, Y)::Ob
        cUR::Hom(munit(), R)
        cUX::Hom(R ⊗ Z, X)
        cUW::Hom(X, W)
        cUY::Hom(R ⊗ W ⊗ Z, Y)
        doX::Hom(munit(), X)
        doZ::Hom(munit(), Z)
        obsX_x::Hom(X, munit())
    end

    prog_overall = @program CausalR_eg3 () begin
        r  = cUR()

        z_for_cx = doZ()
        x_obs = cUX(r, z_for_cx)
        obsX_x(x_obs)

        x_do = doX()
        w    = cUW(x_do)

        z_for_y = doZ()
        y = cUY(r, w, z_for_y)

        return y
    end

    display_var = Set([:Z, :X, :W, :Y])

    expected_lines = [
        "---- Step 4.1 on R-fragment #1 ----",
        "Fragment 1, var Z, case2 = case2b_rewritten",
        "Fragment 1, var W, case2 = skip_no_X_input",
        "Fragment 1, var X, case1 = case1b_detected",
        "Fragment 1, var Y, case2 = skip_no_X_input",
        "---- Step 4.1 on R-fragment #2 ----",
        "Fragment 2, var W, case1 = case1a_do_nothing",
        "Fragment 2, var X, case2 = case2b_rewritten",
    ]

    out = capture_stdout() do
        id_cf_step41!(prog_overall, display_var; verbose=true)
    end

    for line in expected_lines
        @test occursin(line, out)
    end
end

