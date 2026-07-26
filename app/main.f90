program test_graphene
    use types
    use config
    implicit none

    type(AppConfig) :: cfg

    call ReadConfigFromStdin(cfg)

end program test_graphene