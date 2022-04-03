include("elm_line1.jl")

function test_main()
  N = 10
  ratio = 0.1
  xe = elm_line1( 0.0, 1.0, N, ratio )
  for i = 1:N+1
    @printf("%5d %18.10f\n", i, xe[i])
  end
end


test_main()
