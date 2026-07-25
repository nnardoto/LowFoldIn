program test_graphene
    use types
    use toyModels
    implicit none

    type(ToyModel)   :: graphene
    type(FileFormat) :: test_io

    print *, "1. Criando o modelo de Grafeno..."
    call graphene%InitializeGraphene(e0=0.0d0, t=-2.7d0)
    
    print *, "2. Exibindo dados do modelo original:"
    call graphene%DisplayInfo()

    print *, "3. Salvando dados no arquivo binario 'graphene.bin'..."
    call graphene%FileData%WriteFile("graphene.bin")

    print *, "4. Lendo o arquivo 'graphene.bin' para um novo objeto..."
    call test_io%ReadFile("graphene.bin")

    print *, "5. Exibindo dados lidos do arquivo binario:"
    call test_io%DisplayInfo()

    print *, "SUCESSO! O I/O binario e o modelo estao funcionais!"
end program test_graphene