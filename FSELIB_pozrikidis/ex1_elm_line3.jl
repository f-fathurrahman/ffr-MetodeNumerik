include("elm_line3.jl")

function test_main()
  N = 11
  ratio = 1.0
  xe = elm_line3( 0.0, 1.0, N, ratio )
  for i = 1:N+1
    @printf("%5d %18.10f\n", i, xe[i])
  end
end


test_main()
