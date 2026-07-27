program LowFoldIn
    use core
    use types
    use config
    implicit none

    type(AppConfig) :: cfg

    call ReadConfigFromStdin(cfg)
end program LowFoldIn