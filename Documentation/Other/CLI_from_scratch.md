1. npm ```init```ed
2. minimal tsconfig:
    ```json
    {
    "compilerOptions": {
        "baseUrl": ".",
        "target": "ES2017",
        "module": "commonjs",
        "strict": true,
        "esModuleInterop": true,
        "rootDir": ".",
        "outDir": "dist"
    }
    }
3. Watch out for source and destination fields
4. Otherwise just ```tsc``` or ```tsc -w``` in the background and run.
5. Started to reach for ```commander``` for input-parsing.(See ```cli.js```)
6. Create a symlink if needed.